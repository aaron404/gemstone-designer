#version 450
out vec4 out_color;

#define MAX_CUTS 100

uniform vec3    iMouse;
uniform float   iTime;
uniform vec2    iResolution;
uniform int     frame;

uniform float   ior;

uniform int     girdle_facets;
uniform float   girdle_radius;

uniform int     num_cuts;
uniform vec4    cuts[MAX_CUTS];

uniform int     max_bounces;
uniform int     ss;

#define PI 3.14159265
#define MAX_STEPS 20
#define MAX_DIST  20.0
#define SURF_DIST 0.01
#define EPSILON (SURF_DIST * 2.0)

#define deg2rad(x) ((x) * PI / 180.0)


// rotate vector about y axis
vec3 rotate_y(vec3 vec, float angle) {
    float s = sin(angle);
    float c = cos(angle);
    return vec3(vec.x * c - vec.z * s, vec.y, vec.x * s + vec.z * c);
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
    float angle = PI * theta / 180.0;
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

float sdf_gemstone(vec3 p, inout int id) {

    float d;
    float d2;
    id = 0;

    // girdle cut
    if (girdle_facets < 3) {
        d = sdf_cylinder(p, vec3(0.0), girdle_radius);
    } else {
        d = sdf_nprism(p, vec3(0.0), girdle_facets, girdle_radius);
    }

    d = max(sdf_sphere(p, vec3(0.0), 2.5), d);

    float m = (iMouse.y / iResolution.y - 0.5) * 0.5;
    float p2 = (iMouse.x / iResolution.x) * 90.0;

/*
        pub radius: f32,
        pub azimuth: f32,
        pub elevation: f32,
        pub num_facets: f32,
*/
    
    for (int cut=0; cut<num_cuts; cut++) {
        float theta = deg2rad(90.0 - cuts[cut].z);
        d2 = sdf_ncone2(p,
            vec3(0.0, (cuts[cut].x) / cos(theta), 0.0), // radius
            int(cuts[cut].w),
            2.0,
            cuts[cut].z, // cuts[cut].y,
            cuts[cut].y //cuts[cut].z
        );
        //float sub = sdf_nprism(p, vec3(0.0), 12, 2.0 - float(cut) / 10.0);
        if (d2 > d) {
            d = d2;
            id = cut + 1;
        }
        //d = max(sub, d);
    }
    
    return d;
}

float get_dist(vec3 p, inout int id) {
    vec4 s = vec4(3.0, 1.0, 0.0, 1.0);
    
    float dist = MAX_DIST;
    //dist = min(dist, sdf_plane(p, vec3(0.0, 1.0, 0.0), 0.0));
    //dist = min(dist, sdf_sphere(p, s.xyz, s.w));
    //dist = min(dist, sdf_octahedron(p, vec3(-3.0, 1.0, 0.0), 1.0));
    //dist = min(dist, sdf_nprism(p, vec3(0.0, 1.0 + iMouse.y / iResolution.y, 0.0), int(iMouse.x / 25.0), 1.0));
    //dist = min(dist, sdf_ncone2(p, vec3(0.0, 1.0 + iMouse.y / iResolution.y, 0.0), int(iMouse.x / 25.0), -1.0, 30.0, 0.3));
    id = 0;

    float d2;
    d2 = sdf_gemstone(p, id);
    if (d2 < dist) {
        dist = d2;
        id = id;
    } else {
        id = 0;
    }
    //dist = min(dist, sdf_gemstone(p, id));

    return dist;
}

vec3 get_normal(vec3 p) {
    int id;
    float d = get_dist(p, id);
    vec2 e = vec2(0.01, 0.0);
    
    vec3 n = d - vec3(
        get_dist(p - e.xyy, id),
        get_dist(p - e.yxy, id),
        get_dist(p - e.yyx, id));
    
    return normalize(n);
}

float ray_march(vec3 ro, vec3 rd, out int id) {
    // distance from origin
    float dO = 0.0;

    id = 0;
    
    float count = 0.0;
    for (int i=0; i<MAX_STEPS; i++) {
        vec3 p = ro + rd * dO;
        // distance to scene
        float dS = get_dist(p, id);
        dO += dS;
        count += 1.0;
        if (dS < SURF_DIST) {
            break;
        } else if (dO > MAX_DIST) {
            id = -100;
            break;
        }
    }
    
    //return count / 10.0;
    return dO;
}

vec3 skybox(vec3 rd) {
    return mod(rd * 4.0, 1.0);
}

vec3 ray_march3(vec3 ro, vec3 rd) {

    vec3 p = ro;
    int id;

    int bounces = 0;
    for (int i=0; i<MAX_STEPS; i++) {
        float d = get_dist(p, id);
        p += abs(d) * rd;
        if (abs(d) < SURF_DIST) {
            vec3 n = get_normal(p);
            bounces += 1;
            if (d < 0.0) {
                vec3 refr = refract(rd, -n, ior);
                if (length(refr) < 0.01) {
                    // total internal reflection
                    rd = reflect(rd, n);
                    p -= n * EPSILON;
                } else {
                    rd = refr;
                    //p += n * EPSILON;
                    break;
                }
            } else {
                rd = refract(rd, -n, 1.0 / ior);
                return(rd);
                p -= 2.0 * n * EPSILON;
            }
            if (bounces == max_bounces) {
                break;
            }
        }
        if (length(p) > MAX_DIST) {
            break;
        }

    }

    return rd;
}

vec3 gem_march(vec3 ro, vec3 rd) {

    //trace to surface
    vec3 p = ro;
    int id;

    float fresnel = 0.0;
    vec3 refl;

    // march initial hit
    for (int i=0; i<MAX_STEPS; i++) {
        float d = get_dist(p, id);
        p += rd * d;
        if (d < SURF_DIST) {
            // hit
            vec3 n = get_normal(p);
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
        float d = get_dist(p, id);
        p += rd * abs(d);
        if (abs(d) < SURF_DIST) {
            vec3 n = get_normal(p);
            vec3 refr = refract(rd, n, ior);
            if (length(refr) < 0.001) {
                // total internal reflection
                rd = reflect(rd, n);
                p += n * EPSILON;
            } else {
                rd = refr;
                break;
            }
        }
    }


    return mix(rd, refl, fresnel);
    return vec3(fresnel);
}

float get_light(vec3 p) {
    vec3 light_pos = vec3(0.0, 3.0, 0.0);
    light_pos.xz += 6.0 * vec2(sin(iTime), -abs(cos(iTime)));
    vec3 l = normalize(light_pos - p);
    vec3 n = get_normal(p);
    
    float dif = dot(n, l);
    
    int id;
    float d = ray_march(p + n * EPSILON, l, id);
    
    return clamp(min(d, dif), 0.0, 1.0);// / dot(light_pos - p, light_pos - p) * 20.0;
    
}


void main() {
    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = (gl_FragCoord.xy - 0.5 * iResolution.xy) / iResolution.y;
    vec2 ds = 2.0 / iResolution;

    float q = iMouse.x;
    //uv = floor(uv * q) / q;
    
    vec3 col = vec3(0.0);

    float t = iTime * 0.5;
    
    // ray origin (camera pos)
    float r = 6.0;
    vec3 ro = vec3(-r * sin(t), 0.0, -r * cos(t));
    
    // camera looks towards this point
    vec3 lookat = vec3(0.0, 0.0, 0.0);
    
    // ray direction (camera lookdir)
    vec3 rd = normalize(vec3(uv.x, uv.y, 1.0));
    rd = rotate_y(rd, -t);

    int id = 0;
    float d = ray_march(ro, rd, id);

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
            // avg_col += skybox(gem_march(ro + vec3(vec2(i, j) * ds, 0.0), rd));
            avg_col += gem_march(ro + vec3(vec2(i, j) * ds, 0.0), rd);
        }
    }
    out_color = vec4(avg_col / pow(ss, 2.0), 1.0);
    return;


    if (d < MAX_DIST) {
        vec3 refl = reflect(rd, n);
        vec3 refr_r = refract(rd, n, 1.00 / ior);
        vec3 refr_g = refract(rd, n, 1.01 / ior);
        vec3 refr_b = refract(rd, n, 1.02 / ior);
        
        
        out_color = vec4(mod(refr_r * 4.0, 1.0).r, mod(refr_g * 4.0, 1.0).g, mod(refr_b * 4.0, 1.0).b, 1.0);
        return;
    } else {
        out_color = vec4(mod(rd * 4.0, 1.0), 1.0);
        return;
    }


    // Output to screen
    out_color = vec4(col, 1.0);

/*
        pub radius: f32,
        pub azimuth: f32,
        pub elevation: f32,
        pub num_facets: f32,
*/

    if (iMouse.z > 0.0) {
        float f = float(id) / 4.0;
        out_color = vec4(f, 0.0, 0.0, 1.0);
    }
}
/*

float sdf_ncone(vec3 p, vec3 pos, int n, float r) {
    float theta = atan(p.z, p.x) + 0.2;
    float x = (2.0 * PI / float(n));
    theta = abs(mod(theta, x) - x / 2.0);
    float phi   = PI / float(n);
    float mag   = length(p.xz);
    vec2  p2    = mag * vec2(cos(theta), sin(theta));
    
    float t = tan(phi);
    r = r * 1.0 * (pos.y - p.y);
    
    float d = 0.0;
    if (p2.y > r * t) {
        d = length(p2 - vec2(r, r * t));
    } else {
        d = p2.x - r;
    }
    
    return d;
}

vec3 debug_sdf(vec3 p) {
    vec3 p2 = mod(abs(p), 1.0) - 0.025;
    if (p2.x + p2.y < 0.05 || p2.x + p2.z < 0.05 || p2.y + p2.z < 0.05) {
        //return vec3(0.0);
    }
    float d = sdf_gemstone(p);
    float m = mod(abs(d), 1.0);
    if (d < 0.0) {
        return vec3(m, 0, 1.0-m);
    } else {
        return vec3(0, m, 1.0-m);
    }
} */