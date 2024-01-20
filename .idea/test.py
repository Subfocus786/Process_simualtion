import numpy as np

# Get the size of the array from the user
array_size = int(input("Enter the size of the array: "))

# Initialize an empty NumPy array
my_array = np.empty(array_size)

# Loop to get user input for each element
for i in range(array_size):
    my_array[i] = float(input(f"Enter value for element {i}: "))

print("Entered NumPy Array:", my_array)
