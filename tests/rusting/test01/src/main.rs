// This program takes a list of numbers, computes their sum, and prints the average.

fn main() {
    let mut numbers = vec![1, 2, 3, 4, 5];
    let mult_const = 2;
    let avg = compute_avg(&mut numbers, mult_const);
    println!("The average of {:?} is: {}", numbers, avg/mult_const);
}

fn compute_sum(numbers: &mut [i32], mult_const: i32) -> i32 {
    for n in numbers.iter_mut() {
        *n *= mult_const;
    }
    numbers.iter().sum()
}

fn compute_avg(numbers: &mut [i32], mult_const: i32) -> i32 {
    let size = numbers.len() as i32;
    compute_sum(numbers, mult_const) / size
}
