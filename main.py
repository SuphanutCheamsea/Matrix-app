import numpy as np
import tkinter as tk
from tkinter import messagebox

# Class matrix operations
class MatrixOperations:
    def __init__(self, matrix):
        self.matrix = np.array(matrix)
    
    def multiply(self, other):
        if self.matrix.shape[1] == other.matrix.shape[0]:
            return MatrixOperations(np.dot(self.matrix, other.matrix))
        else:
            raise ValueError("Cannot multiply: Incompatible matrix dimensions.")

    def add(self, other):
        if self.matrix.shape == other.matrix.shape:
            return MatrixOperations(self.matrix + other.matrix)
        else:
            raise ValueError("Cannot add: Matrices must have the same dimensions.")

    def subtract(self, other):
        if self.matrix.shape == other.matrix.shape:
            return MatrixOperations(self.matrix - other.matrix)
        else:
            raise ValueError("Cannot subtract: Matrices must have the same dimensions.")

    def multiply_by_scalar(self, scalar):
        return MatrixOperations(self.matrix * scalar)

    def determinant(self):
        if self.matrix.shape[0] == self.matrix.shape[1]:
            return np.linalg.det(self.matrix)
        else:
            raise ValueError("Matrix must be square to compute determinant.")

    def transpose(self):
        return MatrixOperations(self.matrix.T)

    def inverse(self):
        if self.matrix.shape[0] == self.matrix.shape[1]:
            if np.linalg.det(self.matrix) != 0:
                return MatrixOperations(np.linalg.inv(self.matrix))
            else:
                raise ValueError("Matrix is singular and cannot be inverted.")
        else:
            raise ValueError("Matrix must be square to find the inverse.")

    def solve_linear_equations(self, results):
        if self.matrix.shape[0] == self.matrix.shape[1]:
            return np.linalg.solve(self.matrix, results)
        else:
            raise ValueError("Matrix must be square to solve linear equations.")
    
    def to_string(self):
        return str(self.matrix)

# Class GUI
class MatrixApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Matrix Operations")
        self.create_widgets()

    def create_widgets(self):

        tk.Label(self.root, text="Enter Matrix 1:").grid(row=0, column=0, padx=10, pady=5)
        self.matrix_input1 = tk.Text(self.root, height=5, width=30)
        self.matrix_input1.grid(row=0, column=1, padx=10, pady=5)

        tk.Label(self.root, text="Enter Matrix 2 (for operations):").grid(row=1, column=0, padx=10, pady=5)
        self.matrix_input2 = tk.Text(self.root, height=5, width=30)
        self.matrix_input2.grid(row=1, column=1, padx=10, pady=5)

        tk.Label(self.root, text="Enter Result Vector (for linear equations):").grid(row=2, column=0, padx=10, pady=5)
        self.vector_input = tk.Text(self.root, height=3, width=30)
        self.vector_input.grid(row=2, column=1, padx=10, pady=5)

        tk.Label(self.root, text="Select Operation:").grid(row=3, column=0, padx=10, pady=5)
        self.operation_var = tk.StringVar(value="Multiply Matrices")
        self.operation_menu = tk.OptionMenu(self.root, self.operation_var, "Multiply Matrices", "Add Matrices", "Subtract Matrices", "Multiply by Scalar", "Determinant (Matrix 1)", "Transpose (Matrix 1)", "Inverse (Matrix 1)", "Solve Linear Equations")
        self.operation_menu.grid(row=3, column=1, padx=10, pady=5)

        tk.Label(self.root, text="Scalar (for multiplication):").grid(row=4, column=0, padx=10, pady=5)
        self.scalar_entry = tk.Entry(self.root)
        self.scalar_entry.grid(row=4, column=1, padx=10, pady=5)

        self.execute_button = tk.Button(self.root, text="Execute", command=self.execute_operation)
        self.execute_button.grid(row=5, column=0, columnspan=2, pady=10)

        tk.Label(self.root, text="Result:").grid(row=6, column=0, padx=10, pady=5)
        self.result_display = tk.Text(self.root, height=5, width=30)
        self.result_display.grid(row=6, column=1, padx=10, pady=5)

    def execute_operation(self):
        matrix1 = self.get_matrix_from_input(self.matrix_input1)
        matrix2 = self.get_matrix_from_input(self.matrix_input2)
        vector = self.get_matrix_from_input(self.vector_input)

        if matrix1 is None:
            return
        matrix_op1 = MatrixOperations(matrix1)

        try:
            operation = self.operation_var.get()
            result = None

            if operation == "Multiply Matrices":
                if matrix2 is None:
                    messagebox.showerror("Error", "Matrix 2 is required for this operation.")
                    return
                matrix_op2 = MatrixOperations(matrix2)
                result = matrix_op1.multiply(matrix_op2)
            elif operation == "Add Matrices":
                if matrix2 is None:
                    messagebox.showerror("Error", "Matrix 2 is required for this operation.")
                    return
                matrix_op2 = MatrixOperations(matrix2)
                result = matrix_op1.add(matrix_op2)
            elif operation == "Subtract Matrices":
                if matrix2 is None:
                    messagebox.showerror("Error", "Matrix 2 is required for this operation.")
                    return
                matrix_op2 = MatrixOperations(matrix2)
                result = matrix_op1.subtract(matrix_op2)
            elif operation == "Multiply by Scalar":
                scalar = float(self.scalar_entry.get())
                result = matrix_op1.multiply_by_scalar(scalar)
            elif operation == "Determinant (Matrix 1)":
                result = matrix_op1.determinant()
            elif operation == "Transpose (Matrix 1)":
                result = matrix_op1.transpose()
            elif operation == "Inverse (Matrix 1)":
                result = matrix_op1.inverse()
            elif operation == "Solve Linear Equations":
                if vector is None:
                    messagebox.showerror("Error", "Result Vector is required for solving equations.")
                    return
                result = matrix_op1.solve_linear_equations(np.array(vector).flatten())

            self.result_display.delete("1.0", tk.END)
            if isinstance(result, MatrixOperations):
                self.result_display.insert(tk.END, result.to_string())
            else:
                self.result_display.insert(tk.END, str(result))
        except ValueError as e:
            messagebox.showerror("Error", str(e))

    def get_matrix_from_input(self, matrix_input):
        try:
            rows = matrix_input.get("1.0", tk.END).strip().split("\n")
            matrix = [list(map(float, row.split())) for row in rows]
            return matrix
        except ValueError:
            messagebox.showerror("Input Error", "Please enter a valid matrix.")
            return None

# Main
root = tk.Tk()
app = MatrixApp(root)
root.mainloop()
