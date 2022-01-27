use egui::Slider;
//Alias the backend to something less mouthful
use egui_glfw_gl as egui_backend;

use std::time::Instant;

use egui_backend::egui::{vec2, Pos2, Rect};
use egui_glfw_gl::glfw::{Context};

const SCREEN_WIDTH: u32 = 1000;
const SCREEN_HEIGHT: u32 = 600;
mod triangle;

struct MousePos {
    x: f64,
    y: f64,
    z: f32,
}

struct Gem {
    cuts: Vec<triangle::Cut>,
    girdle_radius: f32,
    girdle_facets: u8,
    ior: f32,
    max_bounces: u32,
    ss: u8,
    frame: u32,
}

fn main() {

    // let mut cuts: Vec<triangle::Cut> = Vec::new();
    let mut num_cuts = 0;

    let mut gem = Gem {
        cuts: Vec::new(),
        girdle_radius: 2f32,
        girdle_facets: 12,
        ior: 1.333,
        max_bounces: 16,
        ss: 1,
        frame: 0,
    };

    let mut glfw = glfw::init(glfw::FAIL_ON_ERRORS).unwrap();
    glfw.window_hint(glfw::WindowHint::ContextVersion(3, 2));
    glfw.window_hint(glfw::WindowHint::OpenGlProfile(glfw::OpenGlProfileHint::Core));
    glfw.window_hint(glfw::WindowHint::DoubleBuffer(true));
    glfw.window_hint(glfw::WindowHint::Resizable(true));

    let (mut window, events) = glfw.create_window(SCREEN_WIDTH, SCREEN_HEIGHT, "Gems", glfw::WindowMode::Windowed)
        .expect("Failed to create GLFW window.");

    window.set_char_polling(true);
    window.set_cursor_pos_polling(true);
    window.set_key_polling(true);
    window.set_mouse_button_polling(true);
    window.make_current();
    glfw.set_swap_interval(glfw::SwapInterval::Sync(1));

    gl::load_with(|symbol| window.get_proc_address(symbol) as *const _);

    let mut painter = egui_backend::Painter::new(&mut window, SCREEN_WIDTH, SCREEN_HEIGHT);
    let mut egui_ctx = egui::CtxRef::default();

    let (width, height) = window.get_framebuffer_size();
    let native_pixels_per_point = window.get_content_scale().0;

    let mut egui_input_state = egui_backend::EguiInputState::new(egui::RawInput {
        screen_rect: Some(Rect::from_min_size(
            Pos2::new(0f32, 0f32),
            vec2(width as f32, height as f32) / native_pixels_per_point,
        )),
        pixels_per_point: Some(native_pixels_per_point),
        ..Default::default()
    });
    let start_time = Instant::now();

    //We will draw a crisp white triangle using OpenGL.
    let triangle = triangle::Triangle::new();
    println!("{:?}", triangle.handles);

    let mut mouse_pos = MousePos {x: 0f64, y: 0f64, z: 0f32};

    while !window.should_close() {
        egui_input_state.input.time = Some(start_time.elapsed().as_secs_f64());
        egui_ctx.begin_frame(egui_input_state.input.take());

        //In egui 0.10.0 we seem to be losing the value to pixels_per_point,
        //so setting it every frame now.
        //TODO: Investigate if this is the right way.
        egui_input_state.input.pixels_per_point = Some(native_pixels_per_point);

        //An example of how OpenGL can be used to draw custom stuff with egui
        //overlaying it:
        //First clear the background to something nice.
        unsafe {
            // Clear the screen to black
            //gl::ClearColor(0.455, 0.302, 0.663, 1.0);
            //gl::Clear(gl::COLOR_BUFFER_BIT);
        }
        //Then draw our triangle.
        triangle.draw(triangle::UniformValues {
            mouse: [mouse_pos.x as f32, mouse_pos.y as f32, mouse_pos.z as f32],
            time: start_time.elapsed().as_secs_f32(),
            resolution: [SCREEN_WIDTH as f32, SCREEN_HEIGHT as f32],
            num_cuts: num_cuts,
            cuts: &gem.cuts,
            girdle_facets: gem.girdle_facets,
            girdle_radius: gem.girdle_radius,
            ior: gem.ior,
            max_bounces: gem.max_bounces,
            ss: gem.ss,
            frame: gem.frame,
        });

        gem.frame += 1;

        egui::SidePanel::left("Left Panel").show(&egui_ctx, |ui| {
            ui.add(Slider::new(&mut gem.girdle_radius, 0.0..=4.0).text("Girdle Radius"));
            ui.add(Slider::new(&mut gem.girdle_facets, 0..=32).text("Girdle facets"));
            ui.add(Slider::new(&mut gem.ior, 1.0..=3.0).text("IOR"));
            ui.add(Slider::new(&mut gem.max_bounces, 1..=16).text("max bounces"));
            ui.add(Slider::new(&mut gem.ss, 1..=4).text("supersample"));
            
            let mut i=0;
            for cut in gem.cuts.iter_mut() {
                ui.collapsing(format!("cut {i}"), |ui| {
                    ui.separator();
                    ui.add(Slider::new(&mut cut.num_facets, 3.0..=24.0).text("num facets"));
                    ui.add(Slider::new(&mut cut.azimuth, 0.0..=90.0).text("azimuth"));
                    ui.add(Slider::new(&mut cut.elevation, 0.0..=90.0).text("elevation"));
                    ui.add(Slider::new(&mut cut.radius, 0.0..=2.0).text("radius"));
                });
                i += 1;
            }
            if ui.button("Add cut").clicked() {
                gem.cuts.push(
                    triangle::Cut {
                        radius: 2.0,
                        azimuth: 0f32,
                        elevation: 30f32,
                        num_facets: if gem.girdle_facets < 3 {
                            8.0
                        } else {
                            gem.girdle_facets as f32
                        }}
                );
                num_cuts += 1;
            }
        });

        let (egui_output, paint_cmds) = egui_ctx.end_frame();

        //Handle cut, copy text from egui
        if !egui_output.copied_text.is_empty() {
            egui_backend::copy_to_clipboard(&mut egui_input_state, egui_output.copied_text);
        }

        let paint_jobs = egui_ctx.tessellate(paint_cmds);

        //Note: passing a bg_color to paint_jobs will clear any previously drawn stuff.
        //Use this only if egui is being used for all drawing and you aren't mixing your own Open GL
        //drawing calls with it.
        //Since we are custom drawing an OpenGL Triangle we don't need egui to clear the background.
        painter.paint_jobs(
            None,
            paint_jobs,
            &egui_ctx.texture(),
            native_pixels_per_point,
        );

        for (_, event) in glfw::flush_messages(&events) {
            match event {
                glfw::WindowEvent::Close => window.set_should_close(true),
                glfw::WindowEvent::CursorPos(x, y) => {
                    mouse_pos = MousePos {x, y, z: mouse_pos.z};
                    egui_backend::handle_event(event, &mut egui_input_state);
                }
                glfw::WindowEvent::MouseButton  (button, y, _) => {
                    match button {
                        glfw::MouseButtonLeft => { 
                            if y == glfw::Action::Press {
                                mouse_pos.z = 1.0;
                            } else if y == glfw::Action::Release {
                                mouse_pos.z = 0.0;
                            }
                        },
                        _ => {},
                    }
                    egui_backend::handle_event(event, &mut egui_input_state);
                }
                _ => { egui_backend::handle_event(event, &mut egui_input_state); }
            }
        }
        window.swap_buffers();
        glfw.poll_events();
    }
}
