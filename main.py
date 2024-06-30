from manim import ThreeDScene, ThreeDAxes, Create, DEGREES

# class SquareToCircle(Scene):
#     def construct(self):
#         circle = Circle()
#         square = Square()
#         square.flip(RIGHT)
#         square.rotate(-3 * TAU / 8)
#         circle.set_fill(PINK, opacity=0.5)

#         self.play(Create(square))
#         self.play(Transform(square, circle))
#         self.play(FadeOut(square))


class MyScene(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes()
        axes.add_coordinates()

        self.play(Create(axes), run_time=5)
        self.move_camera(phi=45 * DEGREES, theta=-45 * DEGREES, run_time=3)
        self.begin_ambient_camera_rotation(rate=15 * DEGREES, about="theta")
        self.wait(5)
        self.stop_ambient_camera_rotation()

        self.move_camera(phi=60 * DEGREES, theta=-45 * DEGREES, run_time=3)
        self.wait(3)

        self.move_camera(zoom=0.8)
        self.wait(5)
