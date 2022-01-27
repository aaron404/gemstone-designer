// Draws a simple white triangle
// based on the example from:
// https://github.com/brendanzab/gl-rs/blob/master/gl/examples/triangle.rs

use egui_glfw_gl::gl;
use egui_glfw_gl::gl::types::*;
use std::ffi::CString;
use std::mem;
use std::ptr;
use std::str;

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
}

pub struct Triangle {
    pub vs: GLuint,
    pub fs: GLuint,
    pub program: GLuint,
    pub vao: GLuint,
    pub vbo: GLuint,
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
        // let fs = compile_shader(FS_SRC, gl::FRAGMENT_SHADER);
        let fs = compile_shader(include_str!("gem.frag"), gl::FRAGMENT_SHADER);
        let program = link_program(vs, fs);

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

            Triangle {
                // Create GLSL shaders
                vs,
                fs,
                program,
                vao,
                vbo,
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

            // Draw a triangle from the 3 vertices
            gl::DrawArrays(gl::TRIANGLE_FAN, 0, 4);
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
