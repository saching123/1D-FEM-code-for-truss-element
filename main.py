import numpy as np


def bar_element():
    n = int(input("Enter the number of elements: "))

    element_nodes = []
    for i in range(n + 1):
        element_nodes.append(i + 1)

    # start & end nodes of each element
    snofel = []
    enofel = []
    for i in range(n):
        a = int(input("Enter the start node of element " + str(i + 1) + " : "))
        b = int(input("Enter the end node of element " + str(i + 1) + " : "))
        snofel.append(a)
        enofel.append(b)

    # Get material properties of each element
    element_length = []
    element_area = []
    element_E = []
    element_constant = []
    for i in range(n):
        l = float(input("Enter the length of element " + str(i + 1) + " in mm: "))
        E = float(input("Enter the Young's modulus of element " + str(i + 1) +
                        " in N/mm2: "))
        A = float(input("Enter the area of cross section of element " + str(i + 1)
                        + " in mm2: "))
        element_length.append(l)
        element_E.append(E)
        element_area.append(A)
        element_constant.append(float((E * A) / l))

    # Get force and displacement matrices
    displ_matrix = np.ones((len(element_nodes), 1))
    a = int(input("Enter the number of fixed nodes: "))

    for i in range(a):
        n_n = int(input("Enter the node number of the fixed node: "))
        displ_matrix[(n_n - 1), 0] = 0

    force_matrix = np.zeros((len(element_nodes), 1))
    load_nodes = int(input("Enter the total number of loaded nodes: "))

    for i in range(load_nodes):
        l_n = int(input("Enter the node number of loading: "))
        f = float(input("Enter the axial force at this node in N: "))
        force_matrix[(l_n - 1), 0] = f

    # Get element stiffness matrices
    elstmat = []  # Contains a list of element stiffness matrices
    for i in range(n):
        mat = element_constant[i] * np.array([[1, -1], [-1, 1]])
        elstmat.append(mat)

    # Constructing global stiffness matrix
    gstmat = []  # creating an empty list of global stiffness matrices

    for i in range(n):
        GSM_temp = np.zeros((n + 1, n + 1))
        m = snofel[i]
        p = enofel[i]
        address = [m, p]  # address of columns and rows of gstmat
        for j in range(elstmat[0].ndim):
            for k in range(elstmat[0].ndim):
                a = address[j] - 1
                b = address[k] - 1
                GSM_temp[a, b] = elstmat[i][j][k]
        gstmat.append(GSM_temp)

    global_stiffness_matrix = np.zeros((n + 1, n + 1))

    # print(gstmat)

    for i in range(len(gstmat)):
        global_stiffness_matrix = global_stiffness_matrix + gstmat[i]

    # Matrix reduction
    rcdlist = []
    for i in range(len(element_nodes)):
        if displ_matrix[i, 0] == 0:
            rcdlist.append(i)

    rr_global_stiffness_matrix = np.delete(global_stiffness_matrix, rcdlist, 0)  # row reduction
    rc_global_stiffness_matrix = np.delete(rr_global_stiffness_matrix, rcdlist, 1)  # column reduction
    rgsm = rc_global_stiffness_matrix  # final reduced global stiffness matrix
    rforce_matrix = np.delete(force_matrix, rcdlist, 0)  # reduced force matrix
    rdispl_matrix = np.delete(displ_matrix, rcdlist, 0)  # reduced displacment matrix

    # Solving
    displ_result = np.matmul(np.linalg.inv(rgsm), rforce_matrix)
    rin = 0

    # Print the nodal displacements
    for i in range(len(element_nodes)):
        if displ_matrix[i, 0] == 1:
            displ_matrix[i, 0] = displ_result[rin, 0]
            rin += 1

    for i in range(len(element_nodes)):
        print("Displacement at node " + str(i + 1) + " = " + str(round(displ_matrix[i, 0], 4)) + " mm ")


bar_element()