#pragma once

#ifndef FILE_READING_AND_WRITING
#define FILE_READING_AND_WRITING

#include "base_includes.hpp"

/**
 * @brief Check that the path exists, if not create it
 * 
 * @param filePath 
 */
void create_folder_path(const std::string& filePath)
{
    // get the path to the folder of the file
    fs::path path(filePath);
    fs::path folderPath = path.parent_path();

    // check folder structure exist
    if (!fs::exists(folderPath)) {
        // if not create it
        if (fs::create_directories(folderPath)) {
            std::cout << "Folder structure created: " << folderPath << std::endl;
        // if failed to create, exit the code
        } else {
            std::cerr << "Failed to create folder structure: " << folderPath << std::endl;
            exit(1);
        }
    }
}

/**
 * @brief Generate file path given the file name
 * 
 * @param file_name name of the file
 * @return std::string 
 */
std::string generate_file_path(const std::vector<std::string> (&folder_names), 
    const std::string& file_name)
{
    // get the path to the executable file
    fs::path exe_path = fs::current_path();

    // construct the path to the folder
    //fs::path base_path = exe_path / "data";
    //fs::path base_path = exe_path.parent_path().parent_path() / "data";  // use this when debugging windows
    fs::path base_path = exe_path / "data";  // use this when debugging windows office linux
    for(const std::string& folder_name : folder_names){
        base_path /= folder_name;
    }

    // build final file path
    std::string file_path = (base_path / (file_name + ".txt")).string();

    // check that the path exists, and if not create it
    create_folder_path(file_path);

    return file_path;
}

/**
 * @brief Write vectors of the format
 * {{x1, y1, z1,...},
 *  {x2, y2, z2,...},
 *  ...}
 * to a file in the format
 * x1 y1 z1 ...
 * x2 y2 z2 ...
 * 
 * @param vecs vector of vectors
 * @param folder_name name of the folder to put file in
 * @param file_name name of the file
 */
template<typename T>
void write_vectors_to_file_rowwise(const std::vector<std::vector<T>> vecs,
    const std::vector<std::string> (&folder_names), 
    const std::string& file_name)
{
    // generate path to file
    std::string file_path = generate_file_path(folder_names, file_name);

    // attempt to open file
    std::ofstream out_data(file_path);
    if (!out_data) {
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }

    // write data to file
    for(size_t i{0}; i < vecs.size(); i++) {
        for(size_t j{0}; j < vecs[i].size(); j++){
            out_data << vecs[i][j] << " ";
        }
        out_data << std::endl;
    }
    out_data.close();
}

/**
 * @brief Write vector to file
 * 
 * @param vec vector to write
 * @param folder_name name of the folder to put file in
 * @param file_name name of the file
 */
template<typename T>
void write_single_vector_to_file(const std::vector<T> vec,
    const std::vector<std::string> (&folder_names), 
    const std::string& file_name)
{
    // generate path to file
    std::string file_path = generate_file_path(folder_names, file_name);

    // attempt to open file
    std::ofstream out_data(file_path);
    if (!out_data) {
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }

    // write data to file
    for(size_t j{0}; j < vec.size(); j++) out_data << vec[j] << std::endl;
    out_data.close();
}


#endif