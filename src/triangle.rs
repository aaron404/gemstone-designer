// Draws a simple white triangle
// based on the example from:
// https://github.com/brendanzab/gl-rs/blob/master/gl/examples/triangle.rs

use egui_glfw_gl::gl;
use egui_glfw_gl::gl::types::*;
use exr::prelude::ReadChannels;
use exr::prelude::ReadLayers;
use image::GenericImageView;
use std::ffi::CString;
use std::mem;
use std::ptr;
use std::str;

use crate::triangle::spectrum::NUM_RGB2SPECT_SAMPLES;

extern crate exr;

mod spectrum;

use macros::get_uniform_location;

const VS_SRC: &'static str = "
#version 440
in vec2 position;

void main() {
    gl_Position = vec4(2.0 * position - 1.0, 0.0, 1.0);
}";

const FS_SRC: &'static str = "
#version 440
out vec4 out_color;

uniform vec3 iResolution;
uniform sampler2D fb_texture;

void main() {
    // out_color = vec4(1.0, 0.0, 0.0, 1.0);
    out_color = vec4(texture(fb_texture, gl_FragCoord.xy / iResolution.xy).rgb, 1.0);
    out_color.xyz /= iResolution.z;
    // out_color.xyz *= vec3(10.0);
}";

static VERTEX_DATA: [GLfloat; 8] = [
    0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0
];

static PI: f32 = 3.14159265;
const TRIG_CACHE_SIZE: usize = 1000;

// gets cast to array of f32, so num facets must be f32
pub struct Cut {
    pub radius: f32,
    pub azimuth: f32,
    pub elevation: f32,
    pub num_facets: f32,
}

#[derive(Debug)]
pub struct UniformHandles {
    mouse: GLint,
    time: GLint,
    resolution: GLint,
    frame: GLint,
    num_cuts: GLint,
    cuts: GLint,
    girdle_facets: GLint,
    girdle_radius: GLint,
    table: GLint,
    culet: GLint,
    ior: GLint,
    dispersion: GLint,
    wavelength: GLint,
    max_bounces: GLint,
    ss: GLint,
    skybox_mode: GLint,
    gamma: GLint,
    exposure: GLint,
    render_mode: GLint,
    debug_bool: GLint,
    debug_float: GLint,
}

pub struct UniformValues<'a> {
    pub mouse: [f32; 3],
    pub time: f32,
    pub resolution: [f32; 2],
    pub frame: u32,
    pub num_cuts: u32,
    pub cuts: &'a Vec<Cut>,
    pub girdle_facets: u8,
    pub girdle_radius: f32,
    pub table: f32,
    pub culet: f32,
    pub ior: f32,
    pub dispersion: f32,
    pub wavelength: f32,
    pub max_bounces: u32,
    pub ss: u8,
    pub skybox_mode: u32,
    pub gamma: f32,
    pub exposure: f32,
    pub render_mode: u32,
    pub debug_bool: bool,
    pub debug_float: f32,
}

pub struct Triangle {
    pub vs: GLuint,
    pub fs: GLuint,
    pub program: GLuint,
    pub tonemap_program: GLuint,
    pub vao: GLuint,
    pub vbo: GLuint,
    pub width: i32,
    pub height: i32,
    pub sky_texture: GLuint,
    pub ping_texture: GLuint,
    pub pong_texture: GLuint,
    pub accum_texture_a: GLuint,
    pub accum_texture_b: GLuint,
    pub rgb_spd_texture: GLuint,
    pub handle_fb_texture: i32,
    pub handle_sky_texture: i32,
    pub handle_rgb_spd_texture: i32,
    pub fb_ping: u32,
    pub fb_pong: u32,
    pub fb_accum: u32,
    pub ping: bool,
    pub handles: UniformHandles,
}

pub fn compile_shader(src: &str, ty: GLenum) -> GLuint {
    let shader;
    unsafe {
        shader = gl::CreateShader(ty);
        // Attempt to compile the shader
        let c_str = CString::new(src.as_bytes()).unwrap();
        gl::ShaderSource(shader, 1, &c_str.as_ptr(), ptr::null());
        gl::CompileShader(shader);

        // Get the compile status
        let mut status = gl::FALSE as GLint;
        gl::GetShaderiv(shader, gl::COMPILE_STATUS, &mut status);

        // Fail on error
        if status != (gl::TRUE as GLint) {
            let mut len = 0;
            gl::GetShaderiv(shader, gl::INFO_LOG_LENGTH, &mut len);
            let mut buf = Vec::with_capacity(len as usize);
            buf.set_len((len as usize) - 1); // subtract 1 to skip the trailing null character
            gl::GetShaderInfoLog(
                shader,
                len,
                ptr::null_mut(),
                buf.as_mut_ptr() as *mut GLchar,
            );
            panic!(
                "{}",
                str::from_utf8(&buf).expect("ShaderInfoLog not valid utf8")
            );
        }
    }
    shader
}

pub fn link_program(vs: GLuint, fs: GLuint) -> GLuint {
    unsafe {
        let program = gl::CreateProgram();
        gl::AttachShader(program, vs);
        gl::AttachShader(program, fs);
        gl::LinkProgram(program);
        // Get the link status
        let mut status = gl::FALSE as GLint;
        gl::GetProgramiv(program, gl::LINK_STATUS, &mut status);

        // Fail on error
        if status != (gl::TRUE as GLint) {
            let mut len: GLint = 0;
            gl::GetProgramiv(program, gl::INFO_LOG_LENGTH, &mut len);
            let mut buf = Vec::with_capacity(len as usize);
            buf.set_len((len as usize) - 1); // subtract 1 to skip the trailing null character
            gl::GetProgramInfoLog(
                program,
                len,
                ptr::null_mut(),
                buf.as_mut_ptr() as *mut GLchar,
            );
            panic!(
                "{}",
                str::from_utf8(&buf).expect("ProgramInfoLog not valid utf8")
            );
        }
        program
    }
}

impl Triangle {
    pub fn new(width: i32, height: i32) -> Self {
        // Create Vertex Array Object
        let mut vao = 0;
        let mut vbo = 0;
        let vs = compile_shader(VS_SRC, gl::VERTEX_SHADER);
        let fs = compile_shader(include_str!("gem.frag"), gl::FRAGMENT_SHADER);
        let program = link_program(vs, fs);

        let tonemap_program = link_program(
            compile_shader(VS_SRC, gl::VERTEX_SHADER),
            compile_shader(FS_SRC, gl::FRAGMENT_SHADER),
        );

        let mut fb_ping = 0;
        let mut fb_pong = 0;
        let mut fb_accum = 0;

        let image = exr::prelude::read()
            .no_deep_data()
            .largest_resolution_level()
            .all_channels()
            .all_layers()
            .all_attributes()
            // .from_file("res/photo_studio_loft_hall_16k.exr")
            // .from_file("res/studio_small_09_4k.exr")
            .from_file("res/spiaggia_di_mondello_1k.exr")
            .unwrap();

        println!("file read");
        
        let size = image.attributes.display_window.size;
        let mut sky_data = vec![0f32; size.0 * size.1 * 4];
        println!("vec created");

        for channel in image.layer_data[0].channel_data.list.iter() {
            let c = match channel.name.to_string().as_str() {
                "R" => 0,
                "G" => 1,
                "B" => 2,
                _ => 3,
            };
            for (i, val) in channel.sample_data.values_as_f32().enumerate() {
                sky_data[i * 4 + c] = val;
            }
        }
        println!("data copied");

        let mut trig_data = vec![0f32; TRIG_CACHE_SIZE * 4];
        for i in 0..TRIG_CACHE_SIZE {
            let theta: f32 = 2.0 * PI * (i as f32 / TRIG_CACHE_SIZE as f32);
            trig_data[i * 4 + 0] = theta.sin();
            trig_data[i * 4 + 1] = theta.cos();
            trig_data[i * 4 + 2] = theta.tan();
            trig_data[i * 4 + 3] = 0f32;//f32::atan2()
        }

        // RGB Spectral power density
        let mut rgb_spd_data = vec![0f32; spectrum::NUM_RGB2SPECT_SAMPLES * 4];
        for i in 0..NUM_RGB2SPECT_SAMPLES {
            rgb_spd_data[i * 4 + 0] = 1.0;//spectrum::RGBIllum2SpectRed[i];
            rgb_spd_data[i * 4 + 1] = spectrum::RGBIllum2SpectGreen[i];
            rgb_spd_data[i * 4 + 2] = spectrum::RGBIllum2SpectBlue[i];
            rgb_spd_data[i * 4 + 3] = spectrum::RGBIllum2SpectWhite[i];
        }

        let mut sky_texture = 0;
        let mut trig_texture = 1;
        let mut rgb_spd_texture = 2;
        let mut ping_texture = 0;
        let mut pong_texture = 0;
        let mut accum_texture_a = 0;
        let mut accum_texture_b = 0;
        unsafe {
            gl::GenTextures(1, &mut sky_texture);
            gl::BindTexture(gl::TEXTURE_2D, sky_texture);
            gl::TexImage2D(
                gl::TEXTURE_2D,
                0,
                gl::RGBA16F as i32,
                size.0 as i32,
                size.1 as i32,
                0,
                gl::RGBA,
                gl::FLOAT,
                sky_data.as_mut_ptr() as *const std::os::raw::c_void);
            gl::TexParameteri(gl::TEXTURE_2D, gl::TEXTURE_MIN_FILTER, gl::LINEAR as i32);
            gl::TexParameteri(gl::TEXTURE_2D, gl::TEXTURE_MAG_FILTER, gl::LINEAR as i32);
            gl::BindTexture(gl::TEXTURE_2D, 0);

            gl::GenTextures(1, &mut trig_texture);
            gl::BindTexture(gl::TEXTURE_1D, trig_texture);
            gl::TexImage1D(
                gl::TEXTURE_1D,
                0,
                gl::RGBA as i32,
                TRIG_CACHE_SIZE as i32,
                0,
                gl::RGBA,
                gl::FLOAT,
                trig_data.as_mut_ptr() as *const std::os::raw::c_void);
            gl::TexParameteri(gl::TEXTURE_1D, gl::TEXTURE_MIN_FILTER, gl::LINEAR as i32);
            gl::TexParameteri(gl::TEXTURE_1D, gl::TEXTURE_MAG_FILTER, gl::LINEAR as i32);
            gl::BindTexture(gl::TEXTURE_1D, 0);

            gl::GenTextures(1, &mut rgb_spd_texture);
            gl::BindTexture(gl::TEXTURE_1D, rgb_spd_texture);
            gl::TexImage1D(
                gl::TEXTURE_1D,
                0,
                gl::RGBA as i32,
                NUM_RGB2SPECT_SAMPLES as i32,
                0,
                gl::RGBA,
                gl::DOUBLE,
                rgb_spd_data.as_mut_ptr() as *const std::os::raw::c_void);
            gl::TexParameteri(gl::TEXTURE_1D, gl::TEXTURE_MIN_FILTER, gl::LINEAR as i32);
            gl::TexParameteri(gl::TEXTURE_1D, gl::TEXTURE_MAG_FILTER, gl::LINEAR as i32);
            gl::BindTexture(gl::TEXTURE_1D, 0);

            gl::BindFramebuffer(gl::FRAMEBUFFER, 0);
            gl::FramebufferTexture2D(gl::FRAMEBUFFER, gl::COLOR_ATTACHMENT0, gl::TEXTURE_2D, sky_texture, 0);
            gl::FramebufferTexture1D(gl::FRAMEBUFFER, gl::COLOR_ATTACHMENT1, gl::TEXTURE_1D, trig_texture, 0);
            gl::FramebufferTexture1D(gl::FRAMEBUFFER, gl::COLOR_ATTACHMENT2, gl::TEXTURE_1D, rgb_spd_texture, 0);
        
            gl::GenTextures(1, &mut ping_texture);
            gl::BindTexture(gl::TEXTURE_2D, ping_texture);
            gl::TexImage2D(gl::TEXTURE_2D, 0, gl::RGBA32F as i32, width, height, 0, gl::RGBA, gl::FLOAT, ptr::null());
            gl::TexParameteri(gl::TEXTURE_2D, gl::TEXTURE_MIN_FILTER, gl::LINEAR as i32);
            gl::TexParameteri(gl::TEXTURE_2D, gl::TEXTURE_MAG_FILTER, gl::LINEAR as i32);
            gl::BindTexture(gl::TEXTURE_2D, 0);
            
            gl::GenTextures(1, &mut pong_texture);
            gl::BindTexture(gl::TEXTURE_2D, pong_texture);
            gl::TexImage2D(gl::TEXTURE_2D, 0, gl::RGBA32F as i32, width, height, 0, gl::RGBA, gl::FLOAT, ptr::null());
            gl::TexParameteri(gl::TEXTURE_2D, gl::TEXTURE_MIN_FILTER, gl::LINEAR as i32);
            gl::TexParameteri(gl::TEXTURE_2D, gl::TEXTURE_MAG_FILTER, gl::LINEAR as i32);
            gl::BindTexture(gl::TEXTURE_2D, 0);

            gl::GenTextures(1, &mut accum_texture_a);
            gl::BindTexture(gl::TEXTURE_2D, accum_texture_a);
            gl::TexImage2D(gl::TEXTURE_2D, 0, gl::RGBA as i32, width, height, 0, gl::RGBA, gl::FLOAT, ptr::null());
            gl::TexParameteri(gl::TEXTURE_2D, gl::TEXTURE_MIN_FILTER, gl::LINEAR as i32);
            gl::TexParameteri(gl::TEXTURE_2D, gl::TEXTURE_MAG_FILTER, gl::LINEAR as i32);
            gl::BindTexture(gl::TEXTURE_2D, 0);

            gl::GenTextures(1, &mut accum_texture_b);
            gl::BindTexture(gl::TEXTURE_2D, accum_texture_b);
            gl::TexImage2D(gl::TEXTURE_2D, 0, gl::RGBA as i32, width, height, 0, gl::RGBA, gl::FLOAT, ptr::null());
            gl::TexParameteri(gl::TEXTURE_2D, gl::TEXTURE_MIN_FILTER, gl::LINEAR as i32);
            gl::TexParameteri(gl::TEXTURE_2D, gl::TEXTURE_MAG_FILTER, gl::LINEAR as i32);
            gl::BindTexture(gl::TEXTURE_2D, 0);

            gl::GenFramebuffers(1, &mut fb_ping);
            gl::BindFramebuffer(gl::FRAMEBUFFER, fb_ping);
            gl::FramebufferTexture2D(gl::FRAMEBUFFER, gl::COLOR_ATTACHMENT0, gl::TEXTURE_2D, ping_texture, 0);

            gl::GenFramebuffers(1, &mut fb_pong);
            gl::BindFramebuffer(gl::FRAMEBUFFER, fb_pong);
            gl::FramebufferTexture2D(gl::FRAMEBUFFER, gl::COLOR_ATTACHMENT0, gl::TEXTURE_2D, pong_texture, 0);
            
            gl::GenFramebuffers(1, &mut fb_accum);
            gl::BindFramebuffer(gl::FRAMEBUFFER, fb_accum);
            gl::FramebufferTexture2D(gl::FRAMEBUFFER, gl::COLOR_ATTACHMENT0, gl::TEXTURE_2D, accum_texture_a, 0);
            gl::FramebufferTexture2D(gl::FRAMEBUFFER, gl::COLOR_ATTACHMENT1, gl::TEXTURE_2D, accum_texture_b, 0);
            
            gl::BindFramebuffer(gl::FRAMEBUFFER, 0);
        }

        let handle_mouse = CString::new("iMouse").unwrap();
        let handle_time = CString::new("iTime").unwrap();
        let handle_res = CString::new("iResolution").unwrap();
        let handle_num_cuts = CString::new("num_cuts").unwrap();
        let handle_cuts = CString::new("cuts").unwrap();
        let handle_girdle_facets = CString::new("girdle_facets").unwrap();
        let handle_girdle_radius = CString::new("girdle_radius").unwrap();
        let handle_table = CString::new("table").unwrap();
        let handle_culet = CString::new("culet").unwrap();
        let handle_ior = CString::new("ior").unwrap();
        let handle_dispersion = CString::new("dispersion").unwrap();
        let handle_wavelength = CString::new("wavelength").unwrap();
        let handle_max_bounces = CString::new("max_bounces").unwrap();
        let handle_ss = CString::new("ss").unwrap();
        let handle_frame = CString::new("frame").unwrap();
        let handle_skybox_mode = CString::new("skybox_mode").unwrap();
        let handle_gamma = CString::new("gamma").unwrap();
        let handle_exposure = CString::new("exposure").unwrap();
        let handle_render_mode = CString::new("render_mode").unwrap();
        let handle_debug_bool = CString::new("debug_bool").unwrap();
        let handle_debug_float = CString::new("debug_float").unwrap();

        let handle_sky_texture = CString::new("sky_texture").unwrap();
        let handle_fb_texture = CString::new("fb_texture").unwrap();
        let handle_rgb_spd_texture = CString::new("rgb_spd_texture").unwrap();

        unsafe {
            gl::GenVertexArrays(1, &mut vao);
            gl::GenBuffers(1, &mut vbo);

            let handle_mouse = gl::GetUniformLocation(program, handle_mouse.as_ptr());
            let handle_time = gl::GetUniformLocation(program, handle_time.as_ptr());
            let handle_res = gl::GetUniformLocation(program, handle_res.as_ptr());
            let handle_frame = gl::GetUniformLocation(program, handle_frame.as_ptr());
            let handle_num_cuts = gl::GetUniformLocation(program, handle_num_cuts.as_ptr());
            let handle_cuts = gl::GetUniformLocation(program, handle_cuts.as_ptr());
            let handle_girdle_facets = gl::GetUniformLocation(program, handle_girdle_facets.as_ptr());
            let handle_girdle_radius = gl::GetUniformLocation(program, handle_girdle_radius.as_ptr());
            let handle_table = gl::GetUniformLocation(program, handle_table.as_ptr());
            let handle_culet = gl::GetUniformLocation(program, handle_culet.as_ptr());
            let handle_ior = gl::GetUniformLocation(program, handle_ior.as_ptr());
            let handle_dispersion = gl::GetUniformLocation(program, handle_dispersion.as_ptr());
            let handle_wavelength = gl::GetUniformLocation(program, handle_wavelength.as_ptr());
            let handle_max_bounces = gl::GetUniformLocation(program, handle_max_bounces.as_ptr());
            let handle_ss = gl::GetUniformLocation(program, handle_ss.as_ptr());
            let handle_skybox_mode = gl::GetUniformLocation(program, handle_skybox_mode.as_ptr());
            let handle_gamma = gl::GetUniformLocation(program, handle_gamma.as_ptr());
            let handle_exposure = gl::GetUniformLocation(program, handle_exposure.as_ptr());
            let handle_render_mode = gl::GetUniformLocation(program, handle_render_mode.as_ptr());
            let handle_debug_bool = gl::GetUniformLocation(program, handle_debug_bool.as_ptr());
            let handle_debug_float = gl::GetUniformLocation(program, handle_debug_float.as_ptr());
            
            let handle_sky_texture = gl::GetUniformLocation(program, handle_sky_texture.as_ptr());
            let handle_fb_texture = gl::GetUniformLocation(program, handle_fb_texture.as_ptr());
            let handle_rgb_spd_texture = gl::GetUniformLocation(program, handle_rgb_spd_texture.as_ptr());

            println!("{} {} {} {}", sky_texture, ping_texture, pong_texture, rgb_spd_texture);
            Triangle {
                // Create GLSL shaders
                vs,
                fs,
                program,
                tonemap_program,
                vao,
                vbo,
                width, height,
                sky_texture,
                ping_texture,
                pong_texture,
                accum_texture_a,
                accum_texture_b,
                rgb_spd_texture,
                handle_sky_texture,
                handle_fb_texture,
                handle_rgb_spd_texture,
                fb_ping,
                fb_pong,
                fb_accum,
                ping: true,
                handles: UniformHandles {
                    mouse: handle_mouse,
                    time: handle_time,
                    resolution: handle_res,
                    frame: handle_frame,
                    num_cuts: handle_num_cuts,
                    cuts: handle_cuts,
                    girdle_facets: handle_girdle_facets,
                    girdle_radius: handle_girdle_radius,
                    table: handle_table,
                    culet: handle_culet,
                    ior: handle_ior,
                    dispersion: handle_dispersion,
                    wavelength: handle_wavelength,
                    max_bounces: handle_max_bounces,
                    ss: handle_ss,
                    skybox_mode: handle_skybox_mode,
                    gamma: handle_gamma,
                    exposure: handle_exposure,
                    render_mode: handle_render_mode,
                    debug_bool: handle_debug_bool,
                    debug_float: handle_debug_float,
                },
            }
        }
    }

    pub fn clear(&self) {
        unsafe {
            gl::BindFramebuffer(gl::FRAMEBUFFER, self.fb_ping);
            gl::Clear(gl::COLOR_BUFFER_BIT);
            gl::BindFramebuffer(gl::FRAMEBUFFER, self.fb_pong);
            gl::Clear(gl::COLOR_BUFFER_BIT);
        }
    }

    pub fn draw_realtime(&self, uniforms: &UniformValues) {
        unsafe {
            gl::BindFramebuffer(gl::FRAMEBUFFER, 0);
        }    
        self.draw(&uniforms);

        unsafe {
            gl::BindFramebuffer(gl::FRAMEBUFFER, 0);
        }
    }

    pub fn draw_accurate(&mut self, uniforms: &UniformValues) {
        let draw_fb = if self.ping { self.fb_ping } else { self.fb_pong };
        unsafe {
            gl::BindFramebuffer(gl::FRAMEBUFFER, draw_fb);
        }
        self.draw(&uniforms);
        self.tonemap(&uniforms);
        unsafe {
            gl::BindFramebuffer(gl::FRAMEBUFFER, 0);
        }
        self.ping = !self.ping;
    }

    fn draw(&self, uniforms: &UniformValues) {
        unsafe {
            let read_texture = if self.ping { self.pong_texture } else { self.ping_texture };
            gl::BindVertexArray(self.vao);

            // Create a Vertex Buffer Object and copy the vertex data to it
            gl::BindBuffer(gl::ARRAY_BUFFER, self.vbo);
            gl::BufferData(
                gl::ARRAY_BUFFER,
                (VERTEX_DATA.len() * mem::size_of::<GLfloat>()) as GLsizeiptr,
                mem::transmute(&VERTEX_DATA[0]),
                gl::STATIC_DRAW,
            );

            // Use shader program
            gl::UseProgram(self.program);
            let c_out_color = CString::new("out_color").unwrap();
            gl::BindFragDataLocation(self.program, 0, c_out_color.as_ptr());

            // Specify the layout of the vertex data
            let c_position = CString::new("position").unwrap();
            let pos_attr = gl::GetAttribLocation(self.program, c_position.as_ptr());
            gl::EnableVertexAttribArray(pos_attr as GLuint);
            gl::VertexAttribPointer(
                pos_attr as GLuint,
                2,
                gl::FLOAT,
                gl::FALSE as GLboolean,
                0,
                ptr::null(),
            );

            gl::ActiveTexture(gl::TEXTURE0);
            gl::BindTexture(gl::TEXTURE_2D, self.sky_texture);

            gl::ActiveTexture(gl::TEXTURE1);
            gl::BindTexture(gl::TEXTURE_2D, read_texture);

            gl::Uniform3f(self.handles.mouse, uniforms.mouse[0], uniforms.mouse[1], uniforms.mouse[2]);
            // gl::Uniform1f(self.handles.time, if self.ping { uniforms.time } else { 0f32 });
            gl::Uniform1f(self.handles.time, uniforms.time);
            // gl::Uniform1f(self.handles.time, 1f32);
            gl::Uniform2f(self.handles.resolution, uniforms.resolution[0], uniforms.resolution[1]);
            gl::Uniform1i(self.handles.frame, uniforms.frame as i32);
            gl::Uniform1i(self.handles.num_cuts, uniforms.num_cuts as i32);
            gl::Uniform4fv(self.handles.cuts, uniforms.num_cuts as i32, uniforms.cuts.as_ptr()as *const f32);
            gl::Uniform1i(self.handles.girdle_facets, uniforms.girdle_facets as i32);
            gl::Uniform1f(self.handles.girdle_radius, uniforms.girdle_radius);
            gl::Uniform1f(self.handles.table, uniforms.table);
            gl::Uniform1f(self.handles.culet, uniforms.culet);
            gl::Uniform1f(self.handles.ior, uniforms.ior);
            gl::Uniform1f(self.handles.dispersion, uniforms.dispersion);
            gl::Uniform1f(self.handles.wavelength, uniforms.wavelength);
            gl::Uniform1i(self.handles.max_bounces, uniforms.max_bounces as i32);
            gl::Uniform1i(self.handles.ss, uniforms.ss as i32);
            gl::Uniform1ui(self.handles.skybox_mode, uniforms.skybox_mode);
            gl::Uniform1f(self.handles.gamma, uniforms.gamma);
            gl::Uniform1f(self.handles.exposure, uniforms.exposure);
            gl::Uniform1ui(self.handles.render_mode, uniforms.render_mode);
            gl::Uniform1ui(self.handles.debug_bool, uniforms.debug_bool as u32);
            gl::Uniform1f(self.handles.debug_float, uniforms.debug_float);
            gl::Uniform1i(self.handle_sky_texture, 0);
            gl::Uniform1i(self.handle_fb_texture, 1);

            // Draw a triangle from the 3 vertices
            gl::DrawArrays(gl::TRIANGLE_FAN, 0, 4);
            gl::BindTexture(gl::TEXTURE_2D, 0);
            gl::BindTexture(gl::TEXTURE_1D, 0);


        }
    }

    fn tonemap(&self, uniforms: &UniformValues) {
        let write_texture = if self.ping { self.ping_texture } else { self.pong_texture };

        unsafe {
            // tonemap
            gl::BindFramebuffer(gl::FRAMEBUFFER, 0);

            gl::UseProgram(self.tonemap_program);
            let c_out_color = CString::new("out_color").unwrap();
            gl::BindFragDataLocation(self.program, 0, c_out_color.as_ptr());

            // Specify the layout of the vertex data
            let c_position = CString::new("position").unwrap();
            let pos_attr = gl::GetAttribLocation(self.program, c_position.as_ptr());
            gl::EnableVertexAttribArray(pos_attr as GLuint);
            gl::VertexAttribPointer(
                pos_attr as GLuint,
                2,
                gl::FLOAT,
                gl::FALSE as GLboolean,
                0,
                ptr::null(),
            );

            gl::ActiveTexture(gl::TEXTURE0);
            gl::BindTexture(gl::TEXTURE_2D, write_texture);

            gl::Uniform3f(
                gl::GetUniformLocation(
                    self.tonemap_program, 
                    CString::new("iResolution").unwrap().as_ptr()), uniforms.resolution[0], uniforms.resolution[1], uniforms.max_bounces as f32);
            gl::Uniform1i(
                gl::GetUniformLocation(
                    self.tonemap_program, 
                    CString::new("fb_texture").unwrap().as_ptr()), 0);

            gl::DrawArrays(gl::TRIANGLE_FAN, 0, 4);
            gl::BindTexture(gl::TEXTURE_2D, 0);

            gl::ActiveTexture(gl::TEXTURE0);
            gl::BindTexture(gl::TEXTURE_2D, 0);
            gl::ActiveTexture(gl::TEXTURE1);
            gl::BindTexture(gl::TEXTURE_2D, 0);

        }
    }
}

impl Drop for Triangle {
    fn drop(&mut self) {
        unsafe {
            gl::DeleteProgram(self.program);
            gl::DeleteShader(self.fs);
            gl::DeleteShader(self.vs);
            gl::DeleteBuffers(1, &self.vbo);
            gl::DeleteVertexArrays(1, &self.vao);
        }
    }
}
