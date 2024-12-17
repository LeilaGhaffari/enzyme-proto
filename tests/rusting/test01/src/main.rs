fn main() {
    // ---------------------------------------------------------
    // 2D points
    // ---------------------------------------------------------
    let mut point1 = Point2D { x: 1.1, y: 2.2 };
    let mut point2 = Point2D { x: 2.2, y: 4.4};
    point1.scale(2.0);
    point2.scale(1.0);
    let psum = add_points(point1, point2);
    println!(
        "The 1-norm of the sum of two scaled points is: {}",
        psum.norm1()
    );

    // ---------------------------------------------------------
    // 3x3 matrices
    // ---------------------------------------------------------
    let mut m1 = Mat::new();
    let mut m2 = Mat::new();
    let d1 = [1., 2., 3., 4., 5., 6., 7., 8., 9.];
    let d2: [f64; 9] = d1.map(|x| 2. * x);

    // fill
    m1.fill(&d1);
    m2.fill(&d2);
    m1.print();
    m2.print();

    // add
    let madd = add_mats(&m1, &m2);
    madd.print();

    // multiply
    if !can_multiply(&m1, &m2) {
        println!("The matrices cannot be multiplied.");
    } else {
        let mmult = matmult(&m1, &m2);
        mmult.print();
    }
}

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

struct Mat {
    data: [[f64; 3]; 3] // the first '3' is the number of columns
}

impl Mat {
    fn new() -> Self {
        Self {
            data: [[0.; 3]; 3]
        }
    }
    fn fill(&mut self, d: &[f64]) {
        for i in 0..self.data.len() { // loop over rows
            for j in 0..self.data[i].len() { // loop over columns
                self.data[i][j] = d[i * self.data.len() + j];
            }
        }
    }
    fn print(&self) {
        println!("Matrix elements:");
        for row in self.data.iter() { // Loop over rows
            println!("{:?}", row);    // Print the whole row as a slice
        }
    }
}

fn add_mats(m1: &Mat, m2: &Mat) -> Mat {
    let mut m = Mat::new();
    for i in 0..m.data.len() {
        for j in 0..m.data[i].len() {
            m.data[i][j] = m1.data[i][j] + m2.data[i][j];
        }
    }
    m
}

fn can_multiply(m1: &Mat, m2: &Mat) -> bool {
    m1.data[0].len() == m2.data.len()
}

fn matmult(m1: &Mat, m2: &Mat) -> Mat {
    let mut m = Mat::new();
    let n = m.data.len();
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                m.data[i][j] += m1.data[i][k] * m2.data[k][j];
            }
        }
    }
    m
}
