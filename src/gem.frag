#version 450
out vec4 out_color;

#define MAX_CUTS 100

uniform vec3    iMouse;
uniform float   iTime;
uniform vec2    iResolution;
uniform int     frame;

uniform uint    render_mode;
uniform uint    skybox_mode;
uniform float   gamma;
uniform float   exposure;

uniform float   ior;
uniform int     girdle_facets;
uniform float   girdle_radius;
uniform int     num_cuts;
uniform float   table;
uniform float   culet;
uniform vec4    cuts[MAX_CUTS];

uniform int     max_bounces;
uniform int     ss;

// for A/B testing
uniform bool    debug_bool;


uniform sampler2D sky_texture;

#define PI 3.14159265
#define HALF_PI (3.14159265 / 2.0)
#define MAX_STEPS 20
#define MIN_DIST  4.0
#define MAX_DIST  10.0
#define SURF_DIST 0.005
#define EPSILON (SURF_DIST * 2.0)

#define deg2rad(x) ((x) * PI / 180.0)

const vec3 UP = vec3(0.0, 1.0, 0.0);
const vec3 DOWN = vec3(0.0, -1.0, 0.0);

int eval_count = 0;

// rotate vector about y axis
vec3 rotate_y(vec3 vec, float angle) {
    float s = sin(angle);
    float c = cos(angle);
    return vec3(vec.x * c - vec.z * s, vec.y, vec.x * s + vec.z * c);
}

#define INCLUDE_VIOLET 1
vec3 spectrum(float x, float lambda) {
    // map lambda from 0 -> 1 to 2/3 -> 1 to prevent dips in perceived lightness
    // (see graph of L linked above)
    lambda = (2.0 + lambda) / 3.0;
    
    vec3 abc = vec3(x) - vec3(lambda * 0.5, 0.5, 1.0 - lambda * 0.5);
    vec3 rgb = 0.5 + 0.5 * cos(PI * clamp(abs(2.0 * abc / lambda), 0.0, 1.0));

    if (INCLUDE_VIOLET == 1) {
        float f = 2.0;
        rgb.r += (1.0 - lambda) * 0.25 * (1.0 +  cos(PI * clamp(f * abs(2.0 * (x - 1.0 + lambda / (2.0 * f)) / lambda), 0.0, 1.0)));
    }

    return rgb;
}

float sdf_plane(vec3 p, vec3 n, float h) {
    return dot(p, normalize(n)) - h;
}

float sdf_sphere(vec3 p, vec3 pos, float r) {
    return length(p - pos) - r;
}

float sdf_cylinder(vec3 p, vec3 pos, float r) {
    return length(p.xz - pos.xz) - r;
}

float sdf_nprism(vec3 p, vec3 pos, int n, float r) {
    float theta = atan(p.z, p.x) + PI * 0.5;
    float x = (2.0 * PI / float(n));
    theta = abs(mod(theta, x) - x / 2.0);
    float phi   = PI / float(n);
    float mag   = length(p.xz);
    vec2  p2    = mag * vec2(cos(theta), sin(theta));
    
    float t = tan(phi);
    
    if (p2.y > r * t) {
        return length(p2 - vec2(r, r * t));
    }

    return p2.x - r;
}

float sdf_ncone2(vec3 p, vec3 pos, int n, float r, float theta, float phi) {
    float angle = deg2rad(theta);// PI * theta / 180.0;
    vec3 pn = normalize(vec3(1.0 / tan(angle), 1.0, 0.0));
    p -= pos;
    phi = PI * phi / 180.0;
    vec2 p2 = vec2(cos(phi) * p.x - sin(phi) * p.z, sin(phi) * p.x + cos(phi) * p.z);
    p.xz = p2;
    float a = round(0.5 / PI * atan(p.z, p.x) * float(n)) * PI * 2.0 / float(n);
    float c = cos(a), s = sin(a);
    p.xz = vec2(p.x * c + p.z * s, abs(p.z * c - p.x * s));
    float d = dot(p, pn);
    vec3 h = p - d * pn;
    h.z = min(h.z, h.x * tan(PI / float(n))); // slightly incorrect for n < 6
    return sign(d) * length(p + h * min(sign(h.y), 0.0));
}

// optimized for speed
// TODO: optimize trig calculations
float sdf_ncone3(vec3 p, int n, float radius, float azimuth, float elevation) {
    elevation = -deg2rad(elevation);
    azimuth = deg2rad(azimuth);
    float theta = atan(p.z, p.x);
    float magnitude = length(p.xz);
    theta = mod(theta - PI / float(n) - azimuth, 2.0 * PI / float(n)) - PI / float(n);
    vec2 p2 = magnitude * vec2(sin(theta), cos(theta));
    float d = p2.y - p.y * tan(elevation) - radius / cos(elevation);
    return d * cos(elevation);
}

float sdf_gemstone(vec3 p) {

    eval_count += 1;

    float d;
    float d2;

    // girdle cut
    if (girdle_facets < 3) {
        d = sdf_cylinder(p, vec3(0.0), girdle_radius);
    } else {
        d = sdf_nprism(p, vec3(0.0), girdle_facets, girdle_radius);
    }

    // table and culet
    d = max(d, sdf_plane(p, UP, table));
    d = max(d, sdf_plane(p, DOWN, culet));

    d = max(sdf_sphere(p, vec3(0.0), 2.5), d);

    float m = (iMouse.y / iResolution.y - 0.5) * 0.5;
    float p2 = (iMouse.x / iResolution.x) * 90.0;

    // xyzw rad azi elev n
    for (int cut=0; cut<num_cuts; cut++) {
        vec4 c = cuts[cut];
        float theta = deg2rad(90.0 - cuts[cut].z);
        if (debug_bool) {
            d2 = sdf_ncone3(p, int(c.w), c.x, c.y, c.z);
        } else {
            d2 = sdf_ncone2(p,
                vec3(0.0, (cuts[cut].x) / cos(theta), 0.0), // radius
                int(cuts[cut].w),
                2.0,
                cuts[cut].z, // cuts[cut].y,
                cuts[cut].y //cuts[cut].z
            );
        }

        d = max(d, d2);
    }
    
    return d;
}

float get_dist(vec3 p) {
    vec4 s = vec4(3.0, 1.0, 0.0, 1.0);
    
    float dist = MAX_DIST;
    //dist = min(dist, sdf_plane(p, vec3(0.0, 1.0, 0.0), 0.0));
    //dist = min(dist, sdf_sphere(p, s.xyz, s.w));
    //dist = min(dist, sdf_octahedron(p, vec3(-3.0, 1.0, 0.0), 1.0));
    //dist = min(dist, sdf_nprism(p, vec3(0.0, 1.0 + iMouse.y / iResolution.y, 0.0), int(iMouse.x / 25.0), 1.0));
    //dist = min(dist, sdf_ncone2(p, vec3(0.0, 1.0 + iMouse.y / iResolution.y, 0.0), int(iMouse.x / 25.0), -1.0, 30.0, 0.3));

    float d2;
    d2 = sdf_gemstone(p);
    dist = min(dist, d2);
    if (d2 < dist) {
        dist = d2;
    }
    //dist = min(dist, sdf_gemstone(p, id));

    return dist;
}

vec3 get_normal(vec3 p) {
    float d = get_dist(p);
    vec2 e = vec2(0.01, 0.0);
    
    vec3 n = d - vec3(
        get_dist(p - e.xyy),
        get_dist(p - e.yxy),
        get_dist(p - e.yyx));
    
    return normalize(n);
}

float ray_march(vec3 ro, vec3 rd) {
    // distance from origin
    float dO = 0.0;
    
    float count = 0.0;
    for (int i=0; i<MAX_STEPS; i++) {
        vec3 p = ro + rd * dO;
        // distance to scene
        float dS = get_dist(p);
        dO += dS;
        count += 1.0;
        if (dS < SURF_DIST) {
            break;
        } else if (dO > MAX_DIST) {
            break;
        }
    }
    
    //return count / 10.0;
    return dO;
}

float get_light(vec3 p) {
    vec3 light_pos = vec3(0.0, 5.0, 0.0);
    light_pos.xz += 6.0 * vec2(sin(iTime), -abs(cos(iTime)));
    vec3 l = normalize(light_pos - p);
    vec3 n = get_normal(p);
    
    float dif = dot(n, l);

    return dot(light_pos - p, light_pos - p) * 0.005;
    
}

vec3 skybox(vec3 rd) {
    if (skybox_mode == 0) {
        return mod(rd * 4.0, 1.0);
    } else if (skybox_mode == 1) {
        return rd;
    } else if (skybox_mode == 2) {
        return vec3(mod(rd.x + rd.y + rd.z, 1.0));
    } else {
        rd = normalize(rd);
        float theta = atan(rd.x, rd.z);
        float phi   = atan(-rd.y, length(vec2(rd.x, rd.z)));
        vec2 uv = vec2(0.5 * theta / PI + 0.5, phi / PI + 0.5);
        vec3 col = texture(sky_texture, uv).rgb;
        vec3 mapped = 1.0 - exp(-col * exposure);
        mapped = pow(mapped, vec3(1.0 / gamma));
        return mapped;
    }
}

vec3 gem_march(vec3 ro, vec3 rd) {

    //trace to surface
    vec3 p = ro + MIN_DIST * rd;

    float fresnel = 0.0;
    vec3 refl;

    int num_bounces = 0;

    // march initial hit
    for (int i=0; i<MAX_STEPS; i++) {
        float d = get_dist(p);
        p += rd * d;
        if (d < SURF_DIST) {
            num_bounces += 1;
            // hit
            vec3 n = get_normal(p);
            if (render_mode == 1) {
                return vec3(get_light(p));
            }
            fresnel = pow(1.0 + dot(rd, n), 5.0);
            refl = reflect(rd, n);
            rd = refract(rd, n, 1.0 / ior);
            p -= n * EPSILON;
            break;
        }
        if (d > MAX_DIST) {
            break;
        }
    }

    // march internal reflections
    for (int i=0; i<MAX_STEPS; i++) {
        float d = get_dist(p);
        p += rd * abs(d);
        if (abs(d) < SURF_DIST) {
            num_bounces += 1;
            vec3 n = get_normal(p);
            vec3 refr = refract(rd, -n, ior);
            if (length(refr) < 0.001) {
                // total internal reflection
                rd = reflect(rd, -n);
                p += n * EPSILON;
            } else {
                rd = refr;
                break;
            }
        }
    }


    return mix(skybox(rd), skybox(refl), fresnel);
}


void main() {
    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = (gl_FragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;
    vec2 ds = 2.0 / iResolution;

    float q = iMouse.x;
    //uv = floor(uv * q) / q;
    
    vec3 col = vec3(0.0);

    float t = iTime * 0.05;
    
    // ray origin (camera pos)
    float r = 10.0;
    vec3 ro = vec3(-r * sin(t), 0.0, -r * cos(t));
    
    // camera looks towards this point
    vec3 lookat = vec3(0.0, 0.0, 0.0);
    
    // ray direction (camera lookdir)
    vec3 rd = normalize(vec3(uv.x, uv.y, 1.0));
    rd = rotate_y(rd, -t);

    float d = ray_march(ro, rd);

    // out_color = vec4(vec3(d / 10.0), 1.0);
    // return;
    
    // calculate diffuse lighting
    vec3 p = ro + rd * d;
    float dif = get_light(p);
    
    vec3 n = get_normal(p);

    col = vec3(dif);

    // supersample
    vec3 avg_col = vec3(0.0);
    for (int i=0; i<ss; i++) {
        for (int j=0; j<ss; j++) {
            if (length(uv.xy) > 0.30) {
                avg_col += skybox(rd);
                continue;
            }
            // avg_col += skybox(gem_march(ro + vec3(vec2(i, j) * ds, 0.0), rd));
            avg_col += gem_march(ro + vec3(vec2(i, j) * ds, 0.0), rd);
        }
    }
    if (render_mode == 1) {
        float heat = (eval_count / pow(ss, 2.0)) / ((gamma + 6.0) * (exposure + 6.0));
        col = 1.0 - spectrum(heat * 0.8 + 0.1, 0.0);
        out_color = vec4(col, 1.0);
        return;
    }


    out_color = vec4(avg_col / pow(ss, 2.0), 1.0);
    return;
}
