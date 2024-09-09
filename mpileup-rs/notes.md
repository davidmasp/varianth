# Notes

As this is a learning project too, here I am keeping some notes about the
process.

## Notes on BAM/SAM/VCF standards

### Including/Excluding for read flags

The way that fitering and flags work in the sam/bam standard is though
bitwise operations. These come from single bit flags that cover one
feature per flag.

This bitwise flag can be encoded in a numeric form and these are the values
displayed in the samtools mpileup output and also the internals of noodles. 

There is a helpful [website](https://broadinstitute.github.io/picard/explain-flags.html)
where you can translate from bam flags to numeric
values.

In rust, bit-wise operations work via normal AND/OR operations but applied to
non-booean types. (Actually, the boolean AND, OR is also bitwise I think but
the type only encodes the last bit). Thus, if we want to do a bitwise
AND operation for `1001` and `0111` we can perform it such as ([playground](https://gist.github.com/rust-play/4594faad4727e88271371847c39d86e7)):

```rust
fn main() {
    let exclude_flags: u16 = 9;
    let read_flag: u16 =  7;
    let result = exclude_flags & read_flag;
    if result == 0 {
        println!("read is not excluded");
    }
    else {
        println!("{}", result.to_string());
        println!("read is excluded");
    }
}
```

which returns, read is excluded as they share at least of the flags (the last one corresponding to `1`). Changing `7` to `2 or 0010` then converts all the operations in false and causes the read to be included.

More information [here](https://www.tutorialspoint.com/rust/rust_bitwise_operators.htm).


Samtools, has certain include/exclude defaults. The exclude defaults are:

* 2048 - supp
* 1024 - dup
* 512 - qc_fail

So that's why our default value for exclusion is `3584`.

### Copying certain types (Copy trait)


So I was very confused why the functions didn't consume my number
arguments. It turns out that the functions that implement the
"Copy" trait does not need to be referened as it is very
fast to just clone it.

`String` does not implement the `Copy` trait. See [here](https://gist.github.com/rust-play/c477096518554a315421a5d31687795e).

```rust

fn main() {
    
    let number: u16 = 9;
    
    // does implement Copy
    fn add_1(x: u16) -> u16 {
        x + 1
    }
    
    // does implement Copy
    fn move_number(x: u16) -> u16 {
        x.checked_add(1).unwrap()
    }
    
    // does not implement Copy
    fn move_string(x: String) -> () {
        println!("{}", x);
    }
    
    let aa = add_1(number);
    let bb = add_1(number);
    // this works 
    println!("{}", aa);
    println!("{}", bb);
    let new_number = move_number(number);
    let new_number2 = move_number(number);
    println!("{}", new_number);
    println!("{}", new_number2);
    
    let new_str = "dkfsjkf".to_string();
    move_string(new_str);
    // println!("{}", new_str);
    // this does not work because string is moved.

}
```
