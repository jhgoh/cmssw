#include "Geometry/VeryForwardGeometryBuilder/interface/DetGeomDescBuilder.h"

#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DetectorDescription/DDCMS/interface/DDFilteredView.h"
#include "DetectorDescription/DDCMS/interface/DDDetector.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


/*
 * Generic function to build geo (tree of DetGeomDesc) from old DD compact view.
 */
std::unique_ptr<DetGeomDesc> detgeomdescbuilder::buildDetGeomDescFromCompactView(
const DDCompactView& myCompactView) {
// create DDFilteredView (no filter!!)
DDPassAllFilter filter;
DDFilteredView fv(myCompactView, filter);

if (!fv.firstChild()) {
edm::LogError("PPSGeometryESProducer") << "Filtered view is empty. Cannot build.";
}

// Geo info: sentinel node.
auto geoInfoSentinel = std::make_unique<DetGeomDesc>(fv);

// Construct the tree of children geo info (DetGeomDesc).
do {
// Create node, and add it to the geoInfoSentinel's list.
DetGeomDesc* newGD = new DetGeomDesc(fv);
geoInfoSentinel->addComponent(newGD);
  } while (fv.nextSibling());

  fv.parent();

  edm::LogInfo("PPSGeometryESProducer") << "Successfully built geometry, it has "
                                        << (geoInfoSentinel->components()).size() << " DetGeomDesc nodes.";

  return geoInfoSentinel;
}


/*
 * Generic function to build geo (tree of DetGeomDesc) from DD4hep compact view.
 */
std::unique_ptr<DetGeomDesc> detgeomdescbuilder::buildDetGeomDescFromCompactView(
    const cms::DDCompactView& myCompactView) {
  // create DDFilteredView (no filter!!)
  const cms::DDDetector* mySystem = myCompactView.detector();
  const dd4hep::Volume& worldVolume = mySystem->worldVolume();
  cms::DDFilteredView fv(mySystem, worldVolume);
  if (fv.next(0) == false) {
    edm::LogError("PPSGeometryESProducer") << "Filtered view is empty. Cannot build.";
  }

  const cms::DDSpecParRegistry& allSpecParSections = myCompactView.specpars();
  // Geo info: sentinel node.
  auto geoInfoSentinel = std::make_unique<DetGeomDesc>(fv, allSpecParSections);

  // Construct the tree of children geo info (DetGeomDesc).
  do {
    // Create node, and add it to the geoInfoSentinel's list.
    DetGeomDesc* newGD = new DetGeomDesc(fv, allSpecParSections);
    geoInfoSentinel->addComponent(newGD);
  } while (fv.next(0));

  edm::LogInfo("PPSGeometryESProducer") << "Successfully built geometry, it has "
                                        << (geoInfoSentinel->components()).size() << " DetGeomDesc nodes.";

  return geoInfoSentinel;
}
