#version 440

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
uniform float   dispersion;
uniform float   wavelength;
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
uniform float   debug_float;

uniform sampler2D sky_texture;
uniform sampler2D fb_texture;

#define PI 3.14159265
#define HALF_PI (PI * 0.5)
#define TWO_PI (PI * 2.0)
#define MAX_STEPS 80
#define MIN_DIST  4.0
#define MAX_DIST  10.0
#define SURF_DIST 0.0005
#define EPSILON (SURF_DIST * 2.0)

#define deg2rad(x) ((x) * PI / 180.0)

const vec3 UP = vec3(0.0, 1.0, 0.0);
const vec3 DOWN = vec3(0.0, -1.0, 0.0);
const float CIE_Y_integral = 106.856895;

int eval_count = 0;

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

vec3 srgb2xyz(vec3 srgb) {
    const mat3 M = mat3(
        0.4124564, 0.3575761, 0.1804375,
        0.2126729, 0.7151522, 0.0721750,
        0.0193339, 0.1191920, 0.9503041);

    vec3 srgb_p = pow((srgb + 0.055) / 1.055, vec3(2.4));
    return srgb_p * M;
}

vec3 xyz2srgb(vec3 xyz) {
    const mat3 M = mat3(
         3.2404542, -1.5371385, -0.4985314,
        -0.9692660,  1.8760108,  0.0415560,
         0.0556434, -0.2040259,  1.0572252);
    
    vec3 srgb_p = xyz * M;
    return 1.055 * pow(srgb_p, vec3(1.0 / 2.4)) - 0.055;
}

float G(float x, float mu, float s1, float s2) {
    if (x < mu) {
        return exp(-0.5 * pow(x - mu, 2.0) / pow(s1, 2.0));
    } else {
        return exp(-0.5 * pow(x - mu, 2.0) / pow(s2, 2.0));
    }
}

float X(float lambda) {
    return 1.056 * G(lambda, 599.8, 37.9, 31.0) + 0.362 * G(lambda, 442.0, 16.0, 26.7) - 0.065 * G(lambda, 501.1, 20.4, 26.2);
}

float Y(float lambda) {
    return 0.821 * G(lambda, 568.8, 46.9, 40.5) + 0.286 * G(lambda, 530.9, 16.3, 31.1);
}

float Z(float lambda) {
    return 1.217 * G(lambda, 437.0, 11.8, 36.0) + 0.681 * G(lambda, 459.0, 26.0, 13.8);
}

vec3 XYZ(float lambda) {
    return vec3(X(lambda), Y(lambda), Z(lambda));
}

float cauchy(float lambda, float A, float B) {
    return A + B / (lambda * lambda);
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

vec3 snf_ncone(vec3 p, int n, float radius, float azimuth, float elevation) {
    azimuth = deg2rad(azimuth);
    float theta = atan(p.z, p.x);
    if (theta < 0.0) {
        theta = PI + theta;
    }
    float theta_n = 2.0 * PI / float(n);
    float region = floor((theta + 0.5 * theta_n) / theta_n);
    float theta2 = region * theta_n;
    theta = mod(theta - PI / float(n) - azimuth, 2.0 * PI / float(n)) - PI / float(n);
    vec2 p2 = vec2(sin(theta), cos(theta));
    vec3 nrm = vec3(p2.x, 0.0, 0.0);
    vec2 p3 = vec2(cos(theta2) * nrm.x - sin(theta2) * nrm.y, sin(theta2) * nrm.x + cos(theta2) * nrm.y);
    return normalize(vec3(p3.x, 0.0, p3.y));
}

float sdf_nprism(vec3 p, vec3 pos, int n, float r) {
    float theta_n = 2.0 * PI / float(n);
    float theta = atan(p.z, p.x);// + PI;
    if (theta < 0.0) {
        theta += TWO_PI;
    }
    float x = (2.0 * PI / float(n));
    theta = abs(mod(theta, theta_n) - 0.5 * theta_n);
    float phi   = PI / float(n);
    float mag   = length(p.xz);
    vec2  p2    = mag * vec2(cos(theta), sin(theta));
    
    float t = tan(phi);
    
    if (p2.y > r * t) {
        return length(p2 - vec2(r, r * t));
    }

    return p2.x - r;
}

float sdf_nprism2(vec3 p, int n, float radius) {
    float theta_n = TWO_PI / float(n);
    float theta = atan(p.z, p.x);
    if (theta < 0.0) {
        theta += TWO_PI;
    }
    theta = abs(mod(theta + 0.5 * theta_n, theta_n) - theta_n * 0.5);
    float mag = length(p.xz);
    vec2 p2 = mag * vec2(cos(theta), sin(theta));
    float y = radius * tan(0.5 * theta_n);
    if (p2.y > y) {
        return length(vec2(radius, y) - p2);
    } else {
        return p2.x - radius;
    }
}

vec3 snf_nprism(vec3 p, int n, float radius) {
    float theta_n = TWO_PI / float(n);
    float theta = atan(p.z, p.x);
    // theta = abs(theta);
    if (theta < 0.0) {
        theta += TWO_PI;
    }
    float region = floor((theta + 0.5 * theta_n) / theta_n);
    vec2 t = vec2(sin(theta_n * region), cos(theta_n * region));
    vec2 normal = vec2(1.0, 0.0);
    vec3 rotated = vec3(normal.x * t.y - normal.y * t.x, 0.0, normal.x * t.x + normal.y * t.y);
    return rotated;
}

float sdf_gemstone(vec3 p) {

    eval_count += 1;

    float d;
    float d2;

    // girdle cut
    if (girdle_facets < 3) {
        d = sdf_cylinder(p, vec3(0.0), girdle_radius);
    } else {
        // d = sdf_nprism(p, vec3(0.0), girdle_facets, girdle_radius);
        d = sdf_nprism2(p, girdle_facets, girdle_radius);
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
        d2 = sdf_ncone3(p, int(c.w), c.x, c.y, c.z);

        d = max(d, d2);
    }
    
    return d;
}

vec3 snf_gemstone(vec3 p) {
    float d = sdf_nprism2(p, girdle_facets, girdle_radius);
    float d2 = sdf_sphere(p, vec3(0.0), 2.5);

    // return normalize(p); //normalize(p - 2.5 * normalize(p));

    if (abs(d) <= abs(d2)) {
        //return vec3(0.0, 0.0, 1.0);
        return snf_nprism(p, girdle_facets, girdle_radius);
    } else {
        return normalize(p);// - 2.5 * normalize(p));
    }
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
    float d = sdf_gemstone(p);
    vec2 e = vec2(0.01, 0.0);
    
    vec3 n = d - vec3(
        sdf_gemstone(p - e.xyy),
        sdf_gemstone(p - e.yxy),
        sdf_gemstone(p - e.yyx));
    
    return normalize(n);
}

float ray_march(vec3 ro, vec3 rd) {
    // distance from origin
    float dO = 0.0;
    
    float count = 0.0;
    for (int i=0; i<MAX_STEPS; i++) {
        vec3 p = ro + rd * dO;
        // distance to scene
        float dS = sdf_gemstone(p);
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
        rd = normalize(rd);
        float theta = atan(rd.x, rd.z);
        float phi   = atan(-rd.y, length(vec2(rd.x, rd.z)));
        vec2 uv = vec2(0.5 * theta / PI + 0.5, phi / PI + 0.5);
        vec3 srgb = texture(sky_texture, uv).rgb;
        vec3 exposed = 1.0 - exp(-srgb * exposure);
        vec3 xyz = srgb2xyz(srgb);// * XYZ(debug_float);
        if (max_bounces == 1) {
            return vec3(xyz.x);
        } else if (max_bounces == 2) {
            return vec3(xyz.y);
        } else if (max_bounces == 3) {
            return vec3(xyz.z);
        }
        return xyz;
    } else if (skybox_mode == 1) {
        return rd;
    } else if (skybox_mode == 2) {
        return vec3(mod(rd.x + rd.y + rd.z, 1.0));
    } else if (skybox_mode == 3) {
        return mod(rd * 4.0, 1.0);
    } else {
        return mod(round(rd * 100.0), 4.0) / 3.0;
    }
}

vec3 gem_march(vec3 ro, vec3 rd, float lambda) {

    //trace to surface
    vec3 p = ro + MIN_DIST * rd;

    float fresnel = 0.0;
    vec3 refl;

    int num_bounces = 0;

    float ior_c = cauchy(lambda, ior, dispersion * 1000000.0);

    vec3 n;
    // march initial hit
    for (int i=0; i<MAX_STEPS; i++) {
        float d = get_dist(p);
        p += rd * d;
        if (d < SURF_DIST) {
            num_bounces += 1;
            // hit
            if (debug_bool) {
                n = get_normal(p);
            } else {
                n = snf_gemstone(p);
            }
            //return n;
            if (render_mode == 2) {
                return vec3(-dot(rd, n));
            }
            fresnel = pow(1.0 + dot(rd, n), 5.0);
            refl = reflect(rd, n);
            rd = refract(rd, n, 1.0 / ior_c);
            p -= n * EPSILON;
            break;
        }
        if (d > MAX_DIST) {
            break;
        }
    }

    //return n;

    // march internal reflections
    for (int i=0; i<MAX_STEPS; i++) {
        float d = get_dist(p);
        p += rd * abs(d);
        if (abs(d) < SURF_DIST) {
            num_bounces += 1;
            if (debug_bool) {
                n = get_normal(p);
            } else {
                n = snf_gemstone(p);
            }
            vec3 refr = refract(rd, -n, ior_c);
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

    // float xxx = gl_FragCoord.x / iResolution.x;
    // return;

    float q = iMouse.x;
    //uv = floor(uv * q) / q;
    
    vec3 col = vec3(0.0);

    float t = iTime * 0.25;
    
    // ray origin (camera pos)
    float r = 10.0;
    vec3 ro = vec3(-r * sin(t), 4.0 + 0.0 * sin(4.0*t), -r * cos(t));
    
    // camera target
    vec3 lookat = vec3(0.0, 0.0, 0.0);

    // vector from camera to target
    vec3 lookdir = normalize(lookat - ro);

    // vector pointing to the right in camera space
    vec3 cam_right = cross(lookdir, UP);

    // vector pointing up in camera space
    vec3 cam_up = cross(cam_right, lookdir);
    
    // ray direction (camera lookdir)
    vec3 rd = normalize(uv.x * cam_right + uv.y * cam_up + lookdir);

    float min_freq = 400.0, max_freq = 700.0;

    float lambda = (iMouse.y / iResolution.y) * (max_freq - min_freq) + min_freq;

    int num_bins = 1;
    vec3 color = vec3(0.0);
    vec2 m = iMouse.xy / iResolution.xy;


    vec3 xyz_sum = vec3(0.0);
    if (uv.y < -0.25) {
        float lambda = gl_FragCoord.x / iResolution.x;
        out_color = vec4(XYZ(lambda * 300.0 + 400.0), 1.0);
        return;
    }
    if (uv.y < 0.0) {
        float n = 0.0;
        for (float i=400.0; i<700.0; i+=0.5) {
            xyz_sum += XYZ(i);
            n += 1;
        }
        xyz_sum /= n;
        out_color = vec4(xyz_sum + texture(fb_texture, gl_FragCoord.xy / iResolution.xy).rgb, 1.0);
        return;
    }

    lambda = wavelength;
    if (length(uv.xy) > 0.3) {
        color = skybox(rd);// * XYZ(lambda);
    } else {
        color = gem_march(ro, rd, max(lambda, min_freq));// * XYZ(lambda);
    }

    if (wavelength >= min_freq) {
        color *= XYZ(lambda);
    }

    // tonemap
    color = xyz2srgb(color);
    color = 1.0 - exp(-color * exposure);
    color = pow(color, vec3(1.0 / gamma));

    out_color = vec4(color, 1.0);
    out_color = max(out_color, vec4(0.0));
    out_color.rgb += texture(fb_texture, gl_FragCoord.xy / iResolution.xy).rgb;
    
    // if (wavelength >= min_freq) {
    //     out_color.rgb /= vec3(0.10449, 0.10449, 0.10205);
    // }

    return;
    
    if (render_mode == 1) {
        float heat = (eval_count / pow(ss, 2.0)) / ((gamma) * (exposure));
        col = 1.0 - spectrum(heat * 0.8 + 0.1, 0.0);
        out_color = vec4(col + 0.5, 1.0);//texture(0, );
        return;
    }

    // if (uv.x > 0.0) {
    out_color.rgb += texture(fb_texture, gl_FragCoord.xy / iResolution.xy).rgb;
    // }

    // out_color.xyz = out_color.xyz; //texture(fb_texture, gl_FragCoord.xy / iResolution.xy).rgb;
    // out_color.w = 1.0 + texture(fb_texture, gl_FragCoord.xy / iResolution.xy).w * 0.01;

    return;



    // // supersample
    // vec3 avg_col = vec3(0.0);

    // for (int i=0; i<ss; i++) {
    //     for (int j=0; j<ss; j++) {
    //         // avg_col += skybox(rd);
    //         avg_col += gem_march(ro, rd, ior);
    //         // avg_col += gem_march(ro + vec3(vec2(i, j) * ds, 0.0), rd, ior);
    //     }
    // }

    // avg_col = avg_col / pow(ss, 2.0);
    
    // out_color = vec4(avg_col, 1.0);
    // return;
}
