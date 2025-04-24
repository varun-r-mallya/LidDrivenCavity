import numpy as np
import matplotlib.pyplot as plt


def load_performance_data():
    threads = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    elapsed_time = [4547588, 4580443, 3142521, 2290600, 2028613, 2349594, 2276177, 2192905, 1946464, 1484296, 1914334,
                    1735396, 2040608, 1672974, 1699116, 1654355, 1756129, 1337214, 1685929, 1849847, 1995665]
    return np.array(threads), np.array(elapsed_time)


def plot_performance(threads, elapsed_time):
    plt.figure(figsize=(10, 6))
    plt.plot(threads, elapsed_time, 'bo-')
    plt.xlabel('Number of Threads')
    plt.ylabel('Elapsed Time per Step (ns)')
    plt.title('Performance Analysis: Thread Count vs Elapsed Time')
    plt.grid(True)

    # Calculate and annotate speedup
    baseline_time = elapsed_time[0]
    speedup = baseline_time / elapsed_time

    # Add text box with statistics
    stats_text = f'Max Speedup: {speedup.max():.2f}x\n'
    stats_text += f'Optimal Thread Count: {threads[speedup.argmax()]}'
    plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    threads, times = load_performance_data()
    plot_performance(threads, times)
