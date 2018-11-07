/*******************************************************************************
 * include/util/file_util.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <experimental/filesystem>

[[maybe_unused]] static bool isFile(const std::string& file) {
  return std::experimental::filesystem::exists(file)
         && !std::experimental::filesystem::is_directory(file);
}

[[maybe_unused]] static bool isDirectory(const std::string& dir) {
  return std::experimental::filesystem::is_directory(dir);
}

// returns true, if dir already exists or if dir was created
[[maybe_unused]] static bool createDirectory(const std::string& dir) {
  if(!isDirectory(dir)) {
    return std::experimental::filesystem::create_directories(dir);
  }
  else return true;
}

/******************************************************************************/
