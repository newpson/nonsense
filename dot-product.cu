#include <cstddef>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>

#define BLOCK_SIZE 256

__global__ void dot(const float *a, const float *b, float *result, int n)
{
    __shared__ float shared_mem[BLOCK_SIZE];
    int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < n)
        shared_mem[tid] = a[idx] * b[idx];
    __syncthreads();

    for (int s = blockDim.x/2; s > 0; s >>= 1) {
        if (tid < s)
            shared_mem[tid] += shared_mem[tid + s];
        __syncthreads();
    }

    if (tid == 0)
        atomicAdd(result, shared_mem[0]);
}

int main()
{
    const std::size_t n = 500000;

    auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 engine(seed);
    const int width = 100;
    std::uniform_int_distribution random_value(-width, width);

    std::vector<float> a(n);
    for (float &value: a)
        value = random_value(engine) / 100.0;
    std::vector<float> b(n);
    for (float &value: b)
        value = random_value(engine) / 100.0;

    float *d_a, *d_b, *d_result;
    float result = 0.0;
    const std::size_t num_bytes = n * sizeof(float);
    cudaMalloc(&d_a, num_bytes);
    cudaMalloc(&d_b, num_bytes);
    cudaMalloc(&d_result, sizeof(float));
    cudaMemcpy(d_a, (const void *) a.data(), num_bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, (const void *) b.data(), num_bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_result, &result, sizeof(float), cudaMemcpyHostToDevice);

    cudaEvent_t start_gpu, end_gpu;
    cudaEventCreate(&start_gpu);
    cudaEventCreate(&end_gpu);

    // == (n/BLOCK_SIZE + (n%BLOCK_SIZE != 0))
    std::size_t gridSize = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;

    cudaEventRecord(start_gpu);
    dot<<<gridSize, BLOCK_SIZE>>>(d_a, d_b, d_result, n);

    cudaEventRecord(end_gpu);
    cudaEventSynchronize(end_gpu);
    cudaMemcpy(&result, d_result, sizeof(float), cudaMemcpyDeviceToHost);
    float time_gpu = 0;
    cudaEventElapsedTime(&time_gpu, start_gpu, end_gpu);
    cudaEventDestroy(start_gpu);
    cudaEventDestroy(end_gpu);
    std::cerr << "GPU: " << "" << time_gpu << "ms, result = " << result << std::endl;
    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_result);

    result = 0.0;
    auto start_cpu = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < n; ++i)
        result += a[i] * b[i];
    auto end_cpu = std::chrono::high_resolution_clock::now();
    const float time_cpu = std::chrono::duration_cast<std::chrono::microseconds>(end_cpu - start_cpu).count() / 1000.0f;
    std::cerr << "CPU: " << "" << time_cpu << "ms, result = " << result << std::endl;

    return 0;
}
