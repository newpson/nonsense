#include <iostream>
#include <cstddef>
#include <set>
#include <vector>
#include <random>
#include <chrono>

class CSR_Matrix
{
public:
    using size_t = std::size_t;
    using vector_d = std::vector<double>;
    using vector_i = std::vector<size_t>;

private:
    size_t m_num_rows;
    size_t m_num_columns;
    vector_d m_values;
    vector_i m_columns;
    vector_i m_indices;
public:
    struct Triple
    {
        size_t row;
        size_t column;
        double value;

        // TODO create default trash-filling constructor (?)

        Triple(const size_t row = 0, const size_t column = 0,
               const double value = 0.0)
            : row(row), column(column), value(value)
        {}

        class Compare {
        public:
            bool operator()(const Triple &a, const Triple &b) const
            {
                if (a.row < b.row)
                    return true;
                if (a.row == b.row && a.column < b.column)
                    return true;
                return false;
            }
        };
    };

    CSR_Matrix(const size_t num_rows, const size_t num_columns)
        : m_num_rows(num_rows), m_num_columns(num_columns)
    {}

    size_t num_rows() const { return m_num_rows; }
    size_t num_columns() const { return m_num_columns; }
    const vector_d values() const { return m_values; }
    const vector_i columns() const { return m_columns; }
    const vector_i indices()  const { return m_indices; };

    static CSR_Matrix random(const size_t num_rows, const size_t num_columns,
                             const size_t num_nonzero)
    {
        std::mt19937 engine(112233);
        std::uniform_int_distribution<size_t> random_row(0, num_rows-1);
        std::uniform_int_distribution<size_t> random_column(0, num_columns-1);
        const int width = 100;
        std::uniform_int_distribution random_value(-width, width);

        std::set<Triple, Triple::Compare> triples;
        for (size_t i = 0; i < num_nonzero; ++i) {
            triples.insert(Triple(
                random_row(engine),
                random_column(engine),
                static_cast<double>(random_value(engine)) / static_cast<double>(width)
                ));
        }

        CSR_Matrix matrix(num_rows, num_columns);
        matrix.m_indices.push_back(0); // padding index
        size_t row = 0;
        size_t i = 0;
        for (const Triple &triple: triples) {
            while (row < triple.row) {
                matrix.m_indices.push_back(i);
                ++row;
            }
            matrix.m_values.push_back(triple.value);
            matrix.m_columns.push_back(triple.column);
            ++i;
        }
        matrix.m_indices.push_back(i); // padding index

        return matrix;
    }

    friend std::ostream &operator<<(std::ostream &out, const CSR_Matrix &matrix)
    {
        size_t i = 1;
        for (; i < matrix.m_indices.size(); ++i) {
            size_t j = 0;
            for (size_t col = matrix.m_indices[i-1];
                 j < matrix.num_columns() && col < matrix.m_indices[i]; ++j)
            {
                if (j == matrix.m_columns[col]) {
                    out << matrix.m_values[col];
                    ++col;
                } else {
                    out << 0.0;
                }
                out << "\t";
            }
            for (; j < matrix.num_columns(); ++j)
                out << 0.0 << "\t";
            out << std::endl;
        }
        return out;
    }

    vector_d operator*(const vector_d &a) const
    {
        vector_d result(a.size()); // filled with zeros by default

        for (size_t i = 1, j = 0; i < m_indices.size(); ++i, ++j)
            for (size_t k = m_indices[i-1]; k < m_indices[i]; ++k)
                result[j] += m_values[k] * a[m_columns[k]];

        return result;
    }
};

__global__ void multiply(CSR_Matrix::size_t num_values,
                         CSR_Matrix::size_t num_indices,
                         double *values,
                         CSR_Matrix::size_t *columns,
                         CSR_Matrix::size_t *indices,
                         double *a,
                         double *result)
{
    CSR_Matrix::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < num_indices) {
        CSR_Matrix::size_t j = i - 1;
        result[j] = 0.0;
        for (size_t k = indices[i-1]; k < indices[i]; ++k)
            result[j] += values[k] * a[columns[k]];
    }
}

int main(void)
{
    const CSR_Matrix::size_t num_rows = 5000;
    const CSR_Matrix::size_t num_columns = 5000;
    const CSR_Matrix::size_t num_nonzero = 5000000;

    auto matrix = CSR_Matrix::random(num_rows, num_columns, num_nonzero);

    std::mt19937 engine(112233);
    const int width = 100;
    std::uniform_int_distribution random_value(-width, width);
    CSR_Matrix::vector_d a(num_columns);
    for (double &value: a)
        value = random_value(engine);

    double *d_values;
    const std::size_t bytes_values = matrix.values().size() * sizeof(double);
    CSR_Matrix::size_t *d_columns;
    const std::size_t bytes_columns = matrix.columns().size() * sizeof(CSR_Matrix::size_t);
    CSR_Matrix::size_t *d_indices;
    const std::size_t bytes_indices = matrix.indices().size() * sizeof(CSR_Matrix::size_t);
    double *d_a;
    const std::size_t bytes_a = a.size() * sizeof(double);
    double *d_result;
    const std::size_t bytes_result = bytes_a;
    double *result_gpu = new double[a.size()];

    cudaMalloc(&d_values, bytes_values);
    cudaMalloc(&d_columns, bytes_columns);
    cudaMalloc(&d_indices, bytes_indices);
    cudaMalloc(&d_a, bytes_a);
    cudaMalloc(&d_result, bytes_result);
    cudaMemcpy(d_values, (const void *) matrix.values().data(), bytes_values, cudaMemcpyHostToDevice);
    cudaMemcpy(d_columns, (const void *) matrix.columns().data(), bytes_columns, cudaMemcpyHostToDevice);
    cudaMemcpy(d_indices, (const void *) matrix.indices().data(), bytes_indices, cudaMemcpyHostToDevice);
    cudaMemcpy(d_a, (const void *) a.data(), bytes_a, cudaMemcpyHostToDevice);

    cudaEvent_t start_gpu, end_gpu;
    cudaEventCreate(&start_gpu);
    cudaEventCreate(&end_gpu);
    cudaEventRecord(start_gpu);
    std::size_t block_size = 1024;
    std::size_t grid_size = (matrix.indices().size() + block_size - 1) / block_size;
    multiply<<<grid_size, block_size>>>(matrix.values().size(), matrix.indices().size(),
                                        d_values, d_columns, d_indices, d_a, d_result);
    cudaEventRecord(end_gpu);
    cudaEventSynchronize(end_gpu);
    float time_gpu = 0;
    cudaEventElapsedTime(&time_gpu, start_gpu, end_gpu);
    cudaEventDestroy(start_gpu);
    cudaEventDestroy(end_gpu);

    cudaMemcpy(result_gpu, d_result, bytes_a, cudaMemcpyDeviceToHost);
    std::cerr << "GPU: " << "" << time_gpu << "ms" << std::endl;
    auto start_cpu = std::chrono::high_resolution_clock::now();
    auto result_cpu = matrix * a;
    auto end_cpu = std::chrono::high_resolution_clock::now();
    const float time_cpu = std::chrono::duration_cast<std::chrono::microseconds>(end_cpu - start_cpu).count() / 1000.0f;
    std::cerr << "CPU: " << time_cpu << "ms" << std::endl;

    cudaFree(d_values);
    cudaFree(d_columns);
    cudaFree(d_indices);
    cudaFree(d_a);
    cudaFree(d_result);
    delete[] result_gpu;
    return 0;
}
