#include <cstdint>
#include <cstddef>
#include "Event.h"

Event::Event(uint8_t * input, size_t size) {
    uint8_t * end = input + size;
    number_of_sensors   = *((int32_t*)input); input += sizeof(int32_t);
    number_of_hits      = *((int32_t*)input); input += sizeof(int32_t);
    sensor_Zs           = (int32_t*)input; input += sizeof(int32_t) * number_of_sensors;
    sensor_hits.starts = (int32_t*)input; input += sizeof(int32_t) * number_of_sensors;
    sensor_hits.nums   = (int32_t*)input; input += sizeof(int32_t) * number_of_sensors;
    hit_IDs             = (uint32_t*)input; input += sizeof(uint32_t) * number_of_hits;
    hits.Xs             = (float*)  input; input += sizeof(float)   * number_of_hits;
    hits.Ys             = (float*)  input; input += sizeof(float)   * number_of_hits;
    hits.Zs             = (float*)  input; input += sizeof(float)   * number_of_hits;

}

Event::~Event() {};

void Event::quicksort (float* a, float* b, float* c, unsigned int* d, int start, int end) {
    if (start < end) {
        const int pivot = divide(a, b, c, d, start, end);
        quicksort(a, b, c, d, start, pivot - 1);
        quicksort(a, b, c, d, pivot + 1, end);
    }
}

int Event::divide (float* a, float* b, float* c, unsigned int* d, int start, int end) {
    int left;
    int right;
    float pivot;

    pivot = a[start];
    left = start;
    right = end;

    while (left < right) {
        while (a[right] > pivot) {
            right--;
        }

        while ((left < right) && (a[left] <= pivot)) {
            left++;
        }

        if (left < right) {
            swap(a[left], a[right]);
            swap(b[left], b[right]);
            swap(c[left], c[right]);
            swap(d[left], d[right]);
        }
    }

    swap(a[right], a[start]);
    swap(b[right], b[start]);
    swap(c[right], c[start]);
    swap(d[right], d[start]);

    return right;
}

template<typename T>
void Event::swap (T& a, T& b) {
    T temp = a;
    a = b;
    b = temp;
}
