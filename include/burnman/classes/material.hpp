/*
  TODO: Copyright Notice!
*/
#ifndef BURNMAN_CLASSES_MATERIAL_HPP_INCLUDED
#define BURNMAN_CLASSES_MATERIAL_HPP_INCLUDED

/**
 * Base class for materials.

  TODO
  Python doc:

    Base class for all materials. The main functionality is unroll() which
    returns a list of objects of type :class:`~burnman.Mineral` and their molar
    fractions. This class is available as ``burnman.Material``.

    The user needs to call set_method() (once in the beginning) and set_state()
    before querying the material with unroll() or
 */ 
class Material {

 public:
  // Virtual destructor may be needed for cleanup in derived classes
  virtual ~Material() = default;

};

#endif // BURNMAN_CLASSES_MATERIAL_HPP_INCLUDED