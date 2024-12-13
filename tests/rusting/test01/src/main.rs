// This program takes a list of numbers, computes their sum, and prints the result.

fn main() {
    let numbers = vec![1, 2, 3, 4, 5];
    let sum = compute_sum(&numbers);
    println!("The sum of {:?} is: {}", numbers, sum);
}

fn compute_sum(numbers: &[i32]) -> i32 {
    numbers.iter().sum()
}
