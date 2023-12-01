from typing import List
import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d, Delaunay, delaunay_plot_2d
import matplotlib.pyplot as plt
from .app import main


def read_points(filename: str) -> List[List[int]]:
    """
    Reads the points from the file "in.txt" and returns them as a list of lists.
    
    @return: The list of points.
    @rtype: List[List[int]]

    Time complexity: O(n), where n is the number of points.
    Space complexity: O(n)
    """
    with open(filename, "r") as file:
        n = int(file.readline())
        points = []
        for _ in range(n):
            x = int(file.readline())
            y = int(file.readline())
            points.append([x, y])

    return points

## Sort the coordinates by x ###
def bubleSort(data):
    tempx , tempy = 0, 0
    for i in range(len(data)):
        swapped = False
        for j in range(0, len(data)-i-1):
            if data[j][0] > data[j+1][0]:
                tempx, tempy = data[j][0], data[j][1]
                data[j] = data[j+1]
                data[j+1] = [tempx, tempy]
                swapped = True
        if swapped == False:
            break
    
    return data
## Get the bicectrices ###
def bisec(centrales):
    # Empty list
    bisectrices = [] 
    # for each central do...
    for i in range(len(centrales)): 

        j = i + 1
        
        if j == len(centrales):
            j = 0
        
        # Get the coordinate of central i and i+1
        xa, ya = centrales[i]
        xb, yb = centrales[j]
        #print(i, j)

        # Get the bisectriz from between the centrals
        bisectriz = ((xa + xb) / 2, (ya + yb) / 2)
        bisectrices.append(bisectriz)

    bisectrices = np.array(bisectrices)
    #print(bisectrices)
    return bisectrices

def delaunay(centrales):
    return Delaunay(centrales).simplices

def circumCenter(centrales, tri):
    
    vertex = []
    A = []
    B = []
    C = []
    Ax = []
    Bx = []
    Cx = []
    Ay = []
    By = []
    Cy = []
    
    for i in range(len(tri)):
        contador = 0
        for j in tri[i]:
            if contador == 0:
                A.append(centrales[j])
            elif contador == 1:
                B.append(centrales[j])
            elif contador == 2:
                C.append(centrales[j])
    
            contador = contador + 1
    
    for i in range(len(A)):
        Ax.append(A[i][0])
        Ay.append(A[i][1])
        Bx.append(B[i][0])
        By.append(B[i][1])
        Cx.append(C[i][0])
        Cy.append(C[i][1])
   
    for i in range(len(A)):
        # Calcular el denominador D
        D = 2 * (Ax[i] * (By[i] - Cy[i]) + Bx[i] * (Cy[i] - Ay[i]) + Cx[i] * (Ay[i] - By[i]))
        # Comprobar si D es cercano a cero para evitar divisiones por cero
        if abs(D) < 1e-10:
            raise ValueError("División por cero o denominador cercano a cero. Puede haber colinealidad en los puntos.")
        
        # Calcular las coordenadas del circuncentro
        Ox = (Ax[i]**2 + Ay[i]**2) * (By[i] - Cy[i]) + (Bx[i]**2 + By[i]**2) * (Cy[i] - Ay[i]) + (Cx[i]**2 + Cy[i]**2) * (Ay[i] - By[i])
        Ox /= D
        Oy = (Ax[i]**2 + Ay[i]**2) * (Cx[i] - Bx[i]) + (Bx[i]**2 + By[i]**2) * (Ax[i] - Cx[i]) + (Cx[i]**2 + Cy[i]**2) * (Bx[i] - Ax[i])
        Oy /= D

        vertex.append([Ox, Oy])

    vertex = np.array(vertex)
    vertex = np.unique(vertex, axis=0)
    #vertex = vertex[::-1]
    #print(vertex)

    return np.array(vertex)    

def edges(bisectrices, vertex):

    for i in range(0, 2):
        j = 0
        m =  (vertex[j][1] - bisectrices[i][1]) / (vertex[j][0] - bisectrices[i][0])
        puntos_x = np.linspace(bisectrices[i][0], vertex[j][0], 10)

        # Calcular los correspondientes valores de y utilizando la ecuación de la línea
        puntos_y = m * (puntos_x - bisectrices[i][0]) + bisectrices[i][1]

        plt.plot(puntos_x, puntos_y, linestyle='--', color='red')
    
    for i in range(2, 4):
        j = 1
        m =  (vertex[j][1] - bisectrices[i][1]) / (vertex[j][0] - bisectrices[i][0])
        puntos_x = np.linspace(bisectrices[i][0], vertex[j][0], 10)

        # Calcular los correspondientes valores de y utilizando la ecuación de la línea
        puntos_y = m * (puntos_x - bisectrices[i][0]) + bisectrices[i][1]

        plt.plot(puntos_x, puntos_y, linestyle='--', color='red')
    
    m =  (vertex[0][1] - vertex[1][1]) / (vertex[0][0] - vertex[1][0])
    puntos_x = np.linspace(vertex[1][0], vertex[0][0], 10)

    # Calcular los correspondientes valores de y utilizando la ecuación de la línea
    puntos_y = m * (puntos_x - vertex[1][0]) + vertex[1][1]

    plt.plot(puntos_x, puntos_y, linestyle='--', color='red')

## Build all the necesary to get the voronoy plot
def builVoronoi(centrales):    
    # Find the right order (Idk if this is actually usefull, if not we could use just centrales)
    convexPoints = bubleSort(centrales)
    # Get the bisectices between points
    bisectrices = bisec(centrales)
    # Get Delaunay to get the vertex
    tri = delaunay(convexPoints)
    vertex = circumCenter(centrales, tri)
    

    #edges(bisectrices, vertex, tri)

    return convexPoints, bisectrices, vertex

## Plot the voronoy chart
def plotAll(convexPoints, bisectrices, vertex):

    # Plot the coordinates of bisectrices
    plt.scatter(bisectrices[:, 0], bisectrices[:, 1], label='Bisectriz', color='gray')
    # Plot the central points
    plt.scatter(convexPoints[:,0], convexPoints[:,1], label='Centrales', color='blue')

    plt.scatter(vertex[:, 0], vertex[:, 1], linestyle='--', label='Vertices', color='orange')

    plt.title('Diagrama de Voronoi')
    plt.xlabel('Coordenada X')
    plt.ylabel('Coordenada Y')
    plt.legend()

    vor = Voronoi(convexPoints)
    fig = voronoi_plot_2d(vor)
    print("Vertices:")
    print(vor.vertices)
    print("Bisectrices:")
    print(bisectrices)

    tri = Delaunay(convexPoints)
    fig = delaunay_plot_2d(tri)
    
    plt.show()


def solve_voronoi(points: List[List[int]]) -> None:
    centrales = []
    ### Initial data ###
    option = input("1.- Data from HandTracking\n2.- Data from file\n ")
    if option == "1":  
        main()
        centrales = np.load('centrales.npy')
    elif option == "2":
        centrales = np.array(points)

    print("Centrales:")
    print(centrales)
    # Get the voronoi data
    convexPoints, bisectrices, vertex = builVoronoi(centrales)
    # Plot all the output
    plotAll(centrales, bisectrices, vertex)
