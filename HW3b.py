# I got help from ChatGPT
import math


# Function to calculate the factorial (used in gamma function when alpha is an integer)
def factorial(n):
    if n == 0:
        return 1
    else:
        return n * factorial(n - 1)


# Gamma function approximation using Stirling's formula for non-integer values
def gamma_stirling(x):
    return math.sqrt(2 * math.pi / x) * (x / math.e) ** x


# Gamma function for positive integers or half-integers
def gamma(x):
    if x.is_integer():
        return factorial(int(x) - 1)
    elif (x * 2).is_integer() and x > 0:
        # Using the reflection formula for half-integer factorials
        # Gamma(n+1/2) = (2n)! / (4^n n!) * sqrt(pi)
        n = int(2 * x - 1) // 2
        return factorial(2 * n) / (4 ** n * factorial(n)) * math.sqrt(math.pi)
    else:
        # For non-integer values, use the Stirling's approximation
        return gamma_stirling(x)


# Simpson's rule for numerical integration
def simpsons_rule(f, a, b, n):
    h = (b - a) / n
    k = 0.0
    x = a + h
    for i in range(1, n // 2 + 1):
        k += 4 * f(x)
        x += 2 * h

    x = a + 2 * h
    for i in range(1, n // 2):
        k += 2 * f(x)
        x += 2 * h
    return (h / 3) * (f(a) + f(b) + k)


# t-distribution function to be integrated
def t_distribution(x, m):
    return (1 + x ** 2 / m) ** (-((m + 1) / 2))


# F(z) calculation using the Simpson's rule for numerical integration
def F_z(z, m):
    # Defining the integrand function
    def integrand(u):
        return t_distribution(u, m)

    # Km calculation
    Km = gamma((m + 1) / 2) / (math.sqrt(m * math.pi) * gamma(m / 2))

    # Infinite upper bound for the integration is approximated by a large number
    upper_bound = 10 ** 2  # This can be adjusted for accuracy
    integral_result = simpsons_rule(integrand, -upper_bound, z, 1000)  # 1000 subintervals

    return Km * integral_result


# Let's write a main function to prompt the user for inputs and display the result.
def main():
    # Prompting user for degrees of freedom
    m = int(input("Enter the degrees of freedom: "))
    # Prompting user for the z-value
    z = float(input("Enter the z-value: "))

    # Calculating the probability
    probability = F_z(z, m)

    print(f"The probability for z = {z} with {m} degrees of freedom is: {probability}")


# Run the main function if this script is executed
if __name__ == "__main__":
    main()
