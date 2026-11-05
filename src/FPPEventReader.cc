#include "FPPEventReader.hh"

#include <fstream>
#include <sstream>
#include <stdexcept>

FPPEventReader::FPPEventReader(const std::string& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("Cannot open event file: " + path);

  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '#') continue;
    if (line.rfind("event_id", 0) == 0) continue;
    for (char& c : line) if (c == ',') c = ' ';
    std::istringstream row(line);
    FPPEvent event;
    if (!(row >> event.eventId >> event.weight >> event.xCm >> event.yCm >> event.zCm
              >> event.pxGeV >> event.pyGeV >> event.pzGeV)) {
      throw std::runtime_error("Malformed event row: " + line);
    }
    events_.push_back(event);
  }
  if (events_.empty()) throw std::runtime_error("No events found in: " + path);
}

const FPPEvent& FPPEventReader::event(std::size_t index) const {
  if (index >= events_.size()) throw std::out_of_range("Event index exceeds input sample");
  return events_[index];
}
