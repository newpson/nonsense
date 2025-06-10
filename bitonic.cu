#include <cstddef>
#include <iostream>
#include <chrono>
#include <vector>
#include <random>

__device__ void swap(double *a, double *b)
{
    if (a != b)
    {
        const double temp = *a;
        *a = *b;
        *b = temp;
    }
}

__global__ void bitonic_step(double *arr, std::size_t j, std::size_t k) {
    std::size_t i = threadIdx.x + blockIdx.x * blockDim.x;
    std::size_t ixj = i ^ j;

    if (ixj > i) {
        if ((i & k) == 0) {
            if (arr[i] > arr[ixj])
                swap(arr + i, arr + ixj);
        } else {
            if (arr[i] < arr[ixj])
                swap(arr + i, arr + ixj);
        }
    }
}

bool is_sorted(const double *arr, const std::size_t n) {
    for (std::size_t i = 0; i < n - 1; i++)
        if (arr[i] > arr[i + 1])
            return false;
    return true;
}

void sort_cpu(std::vector<double> &arr) {
    for (int k = 2; k <= arr.size(); k *= 2) {
        for (int j = k/2; j > 0; j /= 2) {
            for (int i = 0; i < arr.size(); i++) {
                int ixj = i ^ j;
                if (ixj > i) {
                    if ((i & k) == 0) {
                        if (arr[i] > arr[ixj])
                            std::swap(arr[i], arr[ixj]);
                    } else {
                        if (arr[i] < arr[ixj])
                            std::swap(arr[i], arr[ixj]);
                    }
                }
            }
        }
    }
}

int main() {
    const std::size_t n = 128;

    std::mt19937 engine(112233);
    const int width = 100;
    std::uniform_int_distribution random_value(-width, width);
    std::vector<double> arr(n);
    for (double &value: arr)
        value = random_value(engine);

    double *d_arr;
    const std::size_t bytes_arr = arr.size() * sizeof(double);
    cudaMalloc(&d_arr, bytes_arr);
    cudaMemcpy(d_arr, (const void *) arr.data(), bytes_arr, cudaMemcpyHostToDevice);
    double *result_gpu = new double[arr.size()];
    const std::size_t block_size = 256;
    const std::size_t grid_size = (n + block_size - 1) / block_size;
    cudaEvent_t start_gpu, end_gpu;
    cudaEventCreate(&start_gpu);
    cudaEventCreate(&end_gpu);
    cudaEventRecord(start_gpu);
    for (std::size_t k = 2; k <= n; k *= 2) {
        for (std::size_t j = k/2; j > 0; j /= 2) {
            bitonic_step<<<grid_size, block_size>>>(d_arr, j, k);
        }
    }
    cudaEventRecord(end_gpu);
    cudaEventSynchronize(end_gpu);
    float time_gpu;
    cudaEventElapsedTime(&time_gpu, start_gpu, end_gpu);
    cudaMemcpy(result_gpu, d_arr, bytes_arr, cudaMemcpyDeviceToHost);
    std::cerr << "GPU: " << time_gpu << "ms, sorted = " << is_sorted(result_gpu, n) << std::endl;
    cudaFree(d_arr);

    auto start_cpu = std::chrono::high_resolution_clock::now();
    sort_cpu(arr);
    auto end_cpu = std::chrono::high_resolution_clock::now();
    const float time_cpu = std::chrono::duration_cast<std::chrono::microseconds>(end_cpu - start_cpu).count() / 1000.0f;
    std::cerr << "CPU: " << time_cpu << "ms, sorted = " << is_sorted(arr.data(), n) << std::endl;

    delete[] result_gpu;

    return 0;
}
