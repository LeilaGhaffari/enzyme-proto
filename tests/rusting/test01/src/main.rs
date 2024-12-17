// This program takes a list of numbers, computes their sum, and prints the average.

struct Point2D {
    x: f64,
    y: f64
}
impl Point2D {
    fn scale(&mut self, mult_const: f64) {
        self.x *= mult_const;
        self.y *= mult_const;
    }
    fn norm1(&self) -> f64 {
        self.x + self.y
    }
}

fn add_points(p1: Point2D, p2: Point2D) -> Point2D {
    Point2D {
        x: p1.x + p2.x,
        y: p1.y + p2.y,
    }
}

fn main() {
    let mut point1 = Point2D { x: 1.1, y: 2.2 };
    let mut point2 = Point2D { x: 2.2, y: 4.4};
    point1.scale(2.0);
    point2.scale(1.0);
    let psum = add_points(point1, point2);
    println!(
        "The 1-norm of the sum of two scaled points is: {}",
        psum.norm1()
    );
}

