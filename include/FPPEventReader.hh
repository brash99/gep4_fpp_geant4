#pragma once

#include "FPPEvent.hh"

#include <string>
#include <vector>

class FPPEventReader {
 public:
  explicit FPPEventReader(const std::string& path);
  const FPPEvent& event(std::size_t index) const;
  std::size_t size() const { return events_.size(); }

 private:
  std::vector<FPPEvent> events_;
};
