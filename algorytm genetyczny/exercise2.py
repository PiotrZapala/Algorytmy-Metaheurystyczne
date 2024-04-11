import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':

    x = np.arange(-10.0, 100.0, 0.5)
    print(x)
    plt.figure()
    plt.plot(x, np.exp(-x) * np.log(1 + np.exp(-x)))
    plt.show()
    