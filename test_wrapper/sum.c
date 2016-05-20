int our_function(int num_numbers, int *numbers) {
    int i;
    int sum;
    for (i = 0; i < num_numbers; i++) {
        sum += numbers[i];
    }
    return sum;
}