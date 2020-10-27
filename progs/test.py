def B_coef(l, m, x):
    """
    """

    if type(l) is not int:
        print("The type of 'l' should be int!")
        sys.exit()
    elif type(m) is not int:
        print("The type of 'm' should be int!")
        sys.exit()

    if m == 0:
        return 0

    if l == m:
        if m == 1:
            return 1
        else:
            return (2*m-1) * m / (m-1) * sqrt(1-x*x) * B_coef(m-1, m-1, x)

    elif l - m == 1:
        return (2*m+1) * x * B_coef(m-1, m-1, x)

    else:
        return ((2*l-1)*x*B_coef(l-1, m, x) - (l-1+m)*B_coef(l-2, m, x)) / (l-m)
