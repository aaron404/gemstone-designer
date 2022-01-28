// Draws a simple white triangle
// based on the example from:
// https://github.com/brendanzab/gl-rs/blob/master/gl/examples/triangle.rs

use egui_glfw_gl::gl;
use egui_glfw_gl::gl::types::*;
use exr::prelude::ReadChannels;
use exr::prelude::ReadLayers;
use std::ffi::CString;
use std::mem;
use std::ptr;
use std::str;

extern crate exr;

const VS_SRC: &'static str = "
#version 450
in vec2 position;

void main() {
    gl_Position = vec4(2.0 * position - 1.0, 0.0, 1.0);
}";

static VERTEX_DATA: [GLfloat; 8] = [
    0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0
];

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
    ior: GLint,
    max_bounces: GLint,
    ss: GLint,
    skybox_mode: GLint,
    gamma: GLint,
    exposure: GLint,
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
    pub ior: f32,
    pub max_bounces: u32,
    pub ss: u8,
    pub skybox_mode: u32,
    pub gamma: f32,
    pub exposure: f32,
}

pub struct Triangle {
    pub vs: GLuint,
    pub fs: GLuint,
    pub program: GLuint,
    pub vao: GLuint,
    pub vbo: GLuint,
    pub sky_texture: GLuint,
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
    pub fn new() -> Self {
        // Create Vertex Array Object
        let mut vao = 0;
        let mut vbo = 0;
        let vs = compile_shader(VS_SRC, gl::VERTEX_SHADER);
        let fs = compile_shader(include_str!("gem.frag"), gl::FRAGMENT_SHADER);
        let program = link_program(vs, fs);

        let image = exr::prelude::read()
            .no_deep_data()
            .largest_resolution_level()
            .all_channels()
            .all_layers()
            .all_attributes()
            .from_file("res/spiaggia_di_mondello_1k.exr")
            .unwrap();
        
        let size = image.attributes.display_window.size;
        let mut data = vec![0f32; size.0 * size.1 * 4];

        for channel in image.layer_data[0].channel_data.list.iter() {
            let c = match channel.name.to_string().as_str() {
                "R" => 0,
                "G" => 1,
                "B" => 2,
                _ => 3,
            };
            for (i, val) in channel.sample_data.values_as_f32().enumerate() {
                data[i * 4 + c] = val;
            }
        }

        let mut sky_texture: GLuint = 1;
        unsafe {
            gl::GenTextures(1, &mut sky_texture);
            gl::BindTexture(gl::TEXTURE_2D, sky_texture);
            gl::TexImage2D(
                gl::TEXTURE_2D,
                0,
                gl::RGBA as i32,
                size.0 as i32,
                size.1 as i32,
                0,
                gl::RGBA,
                gl::FLOAT,
                data.as_mut_ptr() as *const std::os::raw::c_void);

            gl::TexParameteri(gl::TEXTURE_2D, gl::TEXTURE_MIN_FILTER, gl::LINEAR as i32);
            gl::TexParameteri(gl::TEXTURE_2D, gl::TEXTURE_MAG_FILTER, gl::LINEAR as i32);
            gl::BindTexture(gl::TEXTURE_2D, 0);

            gl::BindFramebuffer(gl::FRAMEBUFFER, 0);
            gl::FramebufferTexture2D(gl::FRAMEBUFFER, gl::COLOR_ATTACHMENT0, gl::TEXTURE_2D, sky_texture, 0);
        }

        let handle_mouse = CString::new("iMouse").unwrap();
        let handle_time = CString::new("iTime").unwrap();
        let handle_res = CString::new("iResolution").unwrap();
        let handle_num_cuts = CString::new("num_cuts").unwrap();
        let handle_cuts = CString::new("cuts").unwrap();
        let handle_girdle_facets = CString::new("girdle_facets").unwrap();
        let handle_girdle_radius = CString::new("girdle_radius").unwrap();
        let handle_ior = CString::new("ior").unwrap();
        let handle_max_bounces = CString::new("max_bounces").unwrap();
        let handle_ss = CString::new("ss").unwrap();
        let handle_frame = CString::new("frame").unwrap();
        let handle_skybox_mode = CString::new("skybox_mode").unwrap();
        let handle_gamma = CString::new("gamma").unwrap();
        let handle_exposure = CString::new("exposure").unwrap();

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
            let handle_ior = gl::GetUniformLocation(program, handle_ior.as_ptr());
            let handle_max_bounces = gl::GetUniformLocation(program, handle_max_bounces.as_ptr());
            let handle_ss = gl::GetUniformLocation(program, handle_ss.as_ptr());
            let handle_skybox_mode = gl::GetUniformLocation(program, handle_skybox_mode.as_ptr());
            let handle_gamma = gl::GetUniformLocation(program, handle_gamma.as_ptr());
            let handle_exposure = gl::GetUniformLocation(program, handle_exposure.as_ptr());

            Triangle {
                // Create GLSL shaders
                vs,
                fs,
                program,
                vao,
                vbo,
                sky_texture: sky_texture,
                handles: UniformHandles {
                    mouse: handle_mouse,
                    time: handle_time,
                    resolution: handle_res,
                    frame: handle_frame,
                    num_cuts: handle_num_cuts,
                    cuts: handle_cuts,
                    girdle_facets: handle_girdle_facets,
                    girdle_radius: handle_girdle_radius,
                    ior: handle_ior,
                    max_bounces: handle_max_bounces,
                    ss: handle_ss,
                    skybox_mode: handle_skybox_mode,
                    gamma: handle_gamma,
                    exposure: handle_exposure,
                },
            }
        }
    }
    pub fn draw(&self, uniforms: UniformValues) {
        unsafe {
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

            gl::Uniform3f(self.handles.mouse, uniforms.mouse[0], uniforms.mouse[1], uniforms.mouse[2]);
            gl::Uniform1f(self.handles.time, uniforms.time);
            gl::Uniform2f(self.handles.resolution, uniforms.resolution[0], uniforms.resolution[1]);
            gl::Uniform1i(self.handles.frame, uniforms.frame as i32);
            gl::Uniform1i(self.handles.num_cuts, uniforms.num_cuts as i32);
            gl::Uniform4fv(self.handles.cuts, uniforms.num_cuts as i32, uniforms.cuts.as_ptr()as *const f32);
            gl::Uniform1i(self.handles.girdle_facets, uniforms.girdle_facets as i32);
            gl::Uniform1f(self.handles.girdle_radius, uniforms.girdle_radius);
            gl::Uniform1f(self.handles.ior, uniforms.ior);
            gl::Uniform1i(self.handles.max_bounces, uniforms.max_bounces as i32);
            gl::Uniform1i(self.handles.ss, uniforms.ss as i32);
            gl::Uniform1ui(self.handles.skybox_mode, uniforms.skybox_mode);
            gl::Uniform1f(self.handles.gamma, uniforms.gamma);
            gl::Uniform1f(self.handles.exposure, uniforms.exposure);

            // Draw a triangle from the 3 vertices
            gl::DrawArrays(gl::TRIANGLE_FAN, 0, 4);
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
