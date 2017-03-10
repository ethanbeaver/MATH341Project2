# Code written by Ethan Beaver
import numpy as np
import parser
import logging
import matplotlib.pyplot as plt
# Note: This project is written in Python3. This doesn't guarantee its compatibility with Python2.
#       Please run this project in a Python3 interpreter.

# Additionally: MatPlotLib needs installed through PiP
# Also, tkinter needs to be installed for python3 in order for MatPlotLib to display graphs. It may be installed by
#       default on many OS's.
#       If the graphs aren't showing up, the lack of tkinter is likely the issue.


# Function to perform Romberg integration
#     Takes in input of a function, a left endpoint, a right endpoint, and a tolerance
#     Outputs an approximation for the integral and the h value used
#     Based off pseudocode from the book.
#     Modified to return if tolerance is acceptable.
def romberg(f, a, b, tol):
    n = 100
    h = float(b - a)
    r = np.zeros((2, n+1))
    prev = 0.0
    r[0, 0] = h / 2.0 * (evaluate(f, a) + evaluate(f, b))
    for i in range(2, n+1):
        sm = 0.0
        for k in range(1, int(2 ** (i-2)+1)):
            sm = sm + evaluate(f, a + (k - 0.5) * h)
        r[1, 0] = (1.0 / 2) * (r[0, 0] + h * sm)
        for j in range(1, i):
            r[1, j] = r[1, j - 1] + (r[1, j - 1] - r[0, j - 1]) / (4 ** j - 1)
        # Function now returns if the tolerance is acceptable
        if abs(r[1, 1] - r[0, 0]) < tol and abs(r[0, 0] - prev) < tol and prev != 0:
            return r[1, 1], h
        prev = r[1, 1]
        h /= 2.0
        for j in range(0, i):
            r[0, j] = r[1, j]
    # Function returns integral approximation as well as the
    return r[1, 1], h


# Function to perform Adaptive Quadrature
    # Takes in input of a function, a left endpoint, a right endpoint, a tolerance, and an option to graph or not.
    # Outputs an approximation for the integral and the smallest h value used
    # Based off pseudocode from the book.
    # Modified to keep track of points used as endpoints of Simpson's rule and use my graphing function to plot the
    #       endpoints used by Simpson's method of each subinterval
def adaptive_quad(f, a0, b, tol, graph: bool):
    n = 100
    app = 0.0
    i = 1
    xvals = []
    yvals = []
    n1 = n*n
    # Create arrays to hold data as we perform adaptive quadrature
    tolerance = np.zeros(n1, dtype=np.float)
    a = np.zeros(n1, dtype=np.float)
    h = np.zeros(n1, dtype=np.float)
    f_a = np.zeros(n1, dtype=np.float)
    f_c = np.zeros(n1, dtype=np.float)
    f_b = np.zeros(n1, dtype=np.float)
    s = np.zeros(n1, dtype=np.float)
    l = np.zeros(n1, dtype=np.float)

    # Data for initial interval
    tolerance[i] = 10 * tol
    a[i] = a0
    h[i] = (b - a0)/2
    minh = h[i]
    f_a[i] = evaluate(f, a0)
    f_c[i] = evaluate(f, a0 + h[i])
    f_b[i] = evaluate(f, b)
    s[i] = h[i] * (f_a[i] + 4*f_c[i] + f_b[i])/3
    l[i] = 1
    while i > 0:
        f_d = evaluate(f, a[i] + h[i]/2)
        f_e = evaluate(f, a[i] + 3*h[i]/2)
        s1 = h[i]*(f_a[i] + 4*f_d + f_c[i])/6
        s2 = h[i]*(f_c[i] + 4*f_e + f_b[i])/6
        v1 = a[i]
        v2 = f_a[i]
        v3 = f_c[i]
        v4 = f_b[i]
        v5 = h[i]
        v6 = tolerance[i]
        v7 = s[i]
        v8 = l[i]
        i -= 1
        # If the subinterval is accurate enough, add that data to our approximation of the integral
        if abs(s1 + s2 - v7) < v6:
            app += (s1 + s2)
            # Store the left endpoint of the interval that we used Simpson's rule on for each subinterval
            xvals.append(v1)
            yvals.append(v2)
        else:
            # If we've looped too many times, throw an error
            if v8 >= n:
                logging.error("Level Exceeded")
                return -1
            else:
                # Right half subinterval
                i += 1
                a[i] = v1 + v5
                f_a[i] = v3
                f_c[i] = f_e
                f_b[i] = v4
                h[i] = v5/2
                tolerance[i] = v6/2
                s[i] = s2
                l[i] = v8 + 1
                # Keep track of the smallest h value used
                if h[i] < minh:
                    minh = h[i]

                # Left half subinterval
                i += 1
                a[i] = v1
                f_a[i] = v2
                f_c[i] = f_d
                f_b[i] = v3
                h[i] = h[i - 1]
                tolerance[i] = tolerance[i - 1]
                s[i] = s1
                l[i] = l[i - 1]
                # Keep track of the smallest h value used
                if h[i] < minh:
                    minh = h[i]
    # store the final x and y values (the right endpoint of the original interval)
    xvals.append(b)
    yvals.append(v4)
    # If the user wants to graph the function, graph it
    if graph:
        adaptive_plot(xvals, yvals, f, a0, b)
    return app, minh


# Function to plot data
#     Plots functions based on data from adaptive quadrature
#     Gets inputs of x and y values, as well as the functions and endpoints
#     Graphs the values used in blue and the actual function in a line in red
def adaptive_plot(x, y, f, a, b):
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Adaptive Quadrature')
    x_real = np.linspace(a, b, int(abs(b-a))*50)
    y_real = evaluate(f, x_real)
    plt.plot(x, y, 'bo', x_real, y_real, 'r-')
    plt.show()
    return


# Function to evaluate a user-inputted functions
#     Turns user inputted functions into usable functions in Python
#     Takes input of a user-inputted function and an x value
#     Outputs an evaluation of that function at the x value given
def evaluate(f, x):
    return eval(parser.expr(f).compile())


# Main function
if __name__ == "__main__":
    # Get user input for function, interval, tolerance, and whether or not to graph
    func = input("Please input a math function to use, using x as the variable (python with numpy syntax)\n"
                 "(Ex: x^2*ln(x) would be entered as x**2*np.log(x)): ")
    a = np.float(input("Input left endpoint: "))
    b = np.float(input("Input right endpoint: "))
    tol = 10**np.float(input("Input Tolerance (for 10^-6, enter -6: "))
    try:
        graph = {"yes": True, "no": False}[input('Do you want a graph of Adaptive Quadrature?\n(yes or no): ').lower()]
    except:
        graph = False

    # Compute Romberg integration and Adaptive quadrature and print the results.
    results, h = romberg(func, a, b, tol)
    print("Romberg integration gives an approximation of", results, "with step size", h)
    results, h = adaptive_quad(func, a, b, tol, graph)
    print("Adaptive Quadrature gives an approximation of", results, "with a smallest step size of", h)
