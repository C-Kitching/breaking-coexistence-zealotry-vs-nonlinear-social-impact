#pragma once

#ifndef PROGRESS_BAR
#define PROGRESS_BAR

#include "base_includes.hpp"

// Function to display progress
void display_progress(std::atomic<int>& progress, int total) {
    while (progress < total) {
        int percentage = (100 * progress) / total;
        std::cout << "\rProgress: " << percentage << "%" << std::flush;
        std::this_thread::sleep_for(std::chrono::seconds(1));  // Update every 1s
    }
    std::cout << "\rProgress: 100%" << std::endl;
}

#endif