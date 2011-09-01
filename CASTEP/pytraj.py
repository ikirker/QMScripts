#!/usr/bin/env python
# Gets information from one of my .md files on temperatures, energies, and the Pu-O distances. Probably not useful to anyone else except as an example.
##


import sys
import numpy
import collections
import os


def pds():
  return repr(os.getpid())

#TASTE DEFEAT, PYTHON 2.4
## {{{ http://code.activestate.com/recipes/500261/ (r15)
from operator import itemgetter as _itemgetter
from keyword import iskeyword as _iskeyword
import sys as _sys

def namedtuple(typename, field_names, verbose=False, rename=False):
    """Returns a new subclass of tuple with named fields.

    >>> Point = namedtuple('Point', 'x y')
    >>> Point.__doc__                   # docstring for the new class
    'Point(x, y)'
    >>> p = Point(11, y=22)             # instantiate with positional args or keywords
    >>> p[0] + p[1]                     # indexable like a plain tuple
    33
    >>> x, y = p                        # unpack like a regular tuple
    >>> x, y
    (11, 22)
    >>> p.x + p.y                       # fields also accessable by name
    33
    >>> d = p._asdict()                 # convert to a dictionary
    >>> d['x']
    11
    >>> Point(**d)                      # convert from a dictionary
    Point(x=11, y=22)
    >>> p._replace(x=100)               # _replace() is like str.replace() but targets named fields
    Point(x=100, y=22)

    """

    # Parse and validate the field names.  Validation serves two purposes,
    # generating informative error messages and preventing template injection attacks.
    if isinstance(field_names, basestring):
        field_names = field_names.replace(',', ' ').split() # names separated by whitespace and/or commas
    field_names = tuple(map(str, field_names))
    if rename:
        names = list(field_names)
        seen = set()
        for i, name in enumerate(names):
            if (not min(c.isalnum() or c=='_' for c in name) or _iskeyword(name)
                or not name or name[0].isdigit() or name.startswith('_')
                or name in seen):
                    names[i] = '_%d' % i
            seen.add(name)
        field_names = tuple(names)
    for name in (typename,) + field_names:
        if not min(c.isalnum() or c=='_' for c in name):
            raise ValueError('Type names and field names can only contain alphanumeric characters and underscores: %r' % name)
        if _iskeyword(name):
            raise ValueError('Type names and field names cannot be a keyword: %r' % name)
        if name[0].isdigit():
            raise ValueError('Type names and field names cannot start with a number: %r' % name)
    seen_names = set()
    for name in field_names:
        if name.startswith('_') and not rename:
            raise ValueError('Field names cannot start with an underscore: %r' % name)
        if name in seen_names:
            raise ValueError('Encountered duplicate field name: %r' % name)
        seen_names.add(name)

    # Create and fill-in the class template
    numfields = len(field_names)
    argtxt = repr(field_names).replace("'", "")[1:-1]   # tuple repr without parens or quotes
    reprtxt = ', '.join('%s=%%r' % name for name in field_names)
    template = '''class %(typename)s(tuple):
        '%(typename)s(%(argtxt)s)' \n
        __slots__ = () \n
        _fields = %(field_names)r \n
        def __new__(_cls, %(argtxt)s):
            return _tuple.__new__(_cls, (%(argtxt)s)) \n
        @classmethod
        def _make(cls, iterable, new=tuple.__new__, len=len):
            'Make a new %(typename)s object from a sequence or iterable'
            result = new(cls, iterable)
            if len(result) != %(numfields)d:
                raise TypeError('Expected %(numfields)d arguments, got %%d' %% len(result))
            return result \n
        def __repr__(self):
            return '%(typename)s(%(reprtxt)s)' %% self \n
        def _asdict(self):
            'Return a new dict which maps field names to their values'
            return dict(zip(self._fields, self)) \n
        def _replace(_self, **kwds):
            'Return a new %(typename)s object replacing specified fields with new values'
            result = _self._make(map(kwds.pop, %(field_names)r, _self))
            if kwds:
                raise ValueError('Got unexpected field names: %%r' %% kwds.keys())
            return result \n
        def __getnewargs__(self):
            return tuple(self) \n\n''' % locals()
    for i, name in enumerate(field_names):
        template += '        %s = _property(_itemgetter(%d))\n' % (name, i)
    if verbose:
        print template

    # Execute the template string in a temporary namespace
    namespace = dict(_itemgetter=_itemgetter, __name__='namedtuple_%s' % typename,
                     _property=property, _tuple=tuple)
    try:
        exec template in namespace
    except SyntaxError, e:
        raise SyntaxError(e.message + ':\n' + template)
    result = namespace[typename]

    # For pickling to work, the __module__ variable needs to be set to the frame
    # where the named tuple is created.  Bypass this step in enviroments where
    # sys._getframe is not defined (Jython for example) or sys._getframe is not
    # defined for arguments greater than 0 (IronPython).
    try:
        result.__module__ = _sys._getframe(1).f_globals.get('__name__', '__main__')
    except (AttributeError, ValueError):
        pass

    return result


collections.namedtuple = namedtuple

## end of http://code.activestate.com/recipes/500261/ }}}

#Units:
""" From the NIST website: http://physics.nist.gov/cuu/Constants/index.html """
angstrom_per_bohr    = 0.529177209
femtoseconds_per_aut = 2.4188843265*10**-2
kelvin_per_autemp    = 3.1577464*10**5

def bohr_to_angstrom(val):
  return val*angstrom_per_bohr

def angstrom_to_bohr(val):
  return val/angstrom_per_bohr

def aut_to_fs(val):
  return val*femtoseconds_per_aut
  
def fs_to_aut(val):
  return val/femtoseconds_per_aut

def k_to_autemp(val):
  return val/kelvin_per_autemp
  
def autemp_to_k(val):
  return val*kelvin_per_autemp

def coord_mk(x,y,z):
  return numpy.array([float(x),float(y),float(z)])

class Cell:
  # Okay, I know this is a total cheat, but I am going 
  #  to just assume a cubic cell aligned with the common axes
  #  for now because that's what I need.
  # Details on calculating least image distances in noncubic 
  #  cells can be found here:
  # http://citeseerx.ist.psu.edu/viewdoc/download;jsessionid=D7D50B0BA17F42798C0EF49B36969C52?doi=10.1.1.57.1696&rep=rep1&type=pdf
  # (The Minimum Image Convention in Non-Cubic MD Cells, 
  #   W. Smith, March 29, 1989)
  # 
  # I'm also assuming fixed cell dimensions because argh.
  # (I have put #anothercellassumption everywhere else there are cell assumptions, I hope)
  def __init__(self, a, b, c, angle1, angle2, angle3):
    self.a = a
    self.b = b
    self.c = c
    self.angle1 = angle1
    self.angle2 = angle2
    self.angle3 = angle3
    self.cell_dimensions = coord_mk(a,b,c)
  
  def remap_component(self, distance, cube_side):
    if distance > cube_side/2.0:
      distance = distance - cube_side
    if distance <= -cube_side/2.0:
      distance = distance + cube_side
    return distance
  
  def get_least_image_vector(self, from_vector, to_vector):
    (x,y,z) = numpy.subtract(to_vector, from_vector)
    x = self.remap_component(x, self.a)
    y = self.remap_component(y, self.b)
    z = self.remap_component(z, self.c)
    return coord_mk(x,y,z)
  
  def get_least_image_distance(self, from_vector, to_vector):
    """ Calculates the least image distance between two vectors. """
    result = self.get_least_image_vector(from_vector, to_vector)
    return numpy.linalg.norm(result)

class Trajectory:
  def __init__(self):
    self.atoms = dict()
    self.timesteps = list()
    self.energies = list()
    self.temperatures = list()
    self.cell = None
        
  def number_of_steps(self):
    return len(self.timesteps)
    
class MD_Atom:
  def __init__(self, element, index):
    self.positions = list()
    self.velocities = list()
    self.forces = list()
    self.element = element
    self.index = index
    self.mass = 0
    self.label = element + " " + index
    


def get_all_distances_from_to(trajectory, a_element, b_element=None):
  """Returns data for plotting of all distances of b_element
  away from a_element over time.
  
  If b_element is None, it is set to a_element.
  """
  
  a_filter = lambda x: x.element == a_element
  
  if (b_element == None) or (a_element == b_element):
    b_filter = lambda x: x.element == a_element
    
    list_of_as = filter(a_filter, trajectory.atoms.values())
    list_of_bs = list_of_as
    
    number_of_data = -len(list_of_as) + len(list_of_as) ** 2 
    
  else:
    b_filter = lambda x: x.element == b_element
    
    list_of_as = filter(a_filter, trajectory.atoms.values())
    list_of_bs = filter(b_filter, trajectory.atoms.values())
    number_of_data = len(list_of_as) * len(list_of_bs)
  
  Distance_Set = collections.namedtuple("Distance_Set", 'step time distances')
  distances = list()
  
  for step in range(0, trajectory.number_of_steps() ):
    distances_for_this_step = list()
    for a_atom in list_of_as:
      for b_atom in list_of_bs:
        distances_for_this_step.append( 
          trajectory.cell.get_least_image_distance(
            a_atom.positions[step], b_atom.positions[step]
          )
        )
    distances.append(Distance_Set(step, 
                     trajectory.timesteps[step],
                     distances_for_this_step))
                     
  return distances
  
def graph_distances(distance_data):
  sys.stderr.write("Writing data to distance.data.\n")
  t = open("distance.data"+pds(), "w")
  
  number_of_data = len(distance_data[0].distances)
  fstring = "%f " * (number_of_data)
  for set in distance_data:
    t.write("%f " % aut_to_fs(set.time))
    t.write(fstring % tuple(set.distances))
    #fstring = "%f " * len(set.distances)
    #t.write(fstring % tuple(set.distances))
    t.write("\n")
  
  t.close()
  
  sys.stderr.write("Writing gnuplot input.\n")
  t = open("t.gnuplot"+pds(), "w")
  t.write("set key off; plot ")
  
  for i in range(2,len(distance_data[0].distances)):
    t.write("\"distance.data"+pds()+"\" using 1:%d with lines, " % i)
  t.write("\"distance.data"+pds()+"\" using 1:%d with lines " % (number_of_data+1))
  t.close()
  

def graph_energies(trajectory):
  sys.stderr.write("Writing energies.\n")
  t = open("energies.data"+pds(), "w")
  
  step=0
  t.write("Time TotalE HamiltE KE\n")
  for energy_triple in trajectory.energies:
    t.write("%f %f %f %f\n" % (aut_to_fs(trajectory.timesteps[step]), 
                               energy_triple[0], 
                               energy_triple[1], 
                               energy_triple[2]))
    step=step+1
  t.close()

def graph_temperatures(trajectory):
  sys.stderr.write("Writing temperatures.\n")
  t = open("temperatures.data"+pds(), "w")
  
  step=0
  for temp in trajectory.temperatures:
    t.write("%f %f\n" % (aut_to_fs(trajectory.timesteps[step]), autemp_to_k(temp)))
    step=step+1
  t.close()  

def calculate_rdf(trajectory, from_element, to_element=None):
  # This routine doesn't seem to work correctly
  """Returns a list containing points on the RDF for
  atoms of from_element to atoms of to_element.
  
  For example, "Pu" and "O".
  
  If to_element is None, it is assumed to be the same
  as from_element.
  """
  
  """
  This explanation is for three-dimensional data. To calculate g(r), do the following: 
    * Pick a value of dr 
    * Loop over all values of r that you care about: 
      1) Consider each particle you have in turn. Count all particles that are 
          a distance between r and r + dr away from the particle you're 
          considering. You can think of this as all particles in a spherical 
          shell surrounding the reference particle. The shell has a thickness 
          dr. 
      2) Divide your total count by N, the number of reference particles you 
          considered -- probably the total number of particles in your data. 
      3) Divide this number by 4 pi r^2 dr, the volume of the spherical shell 
          (the surface area 4 pi r^2, multiplied by the small thickness dr). 
          This accounts for the fact that as r gets larger, for trivial reasons 
          you find more particles with the given separation. 
      4) Divide this by the particle number density. This ensures that g(r)=1 
          for data with no structure. In other words, if you just had an 
          arbitrarily placed spherical shell of inner radius r and outer radius 
          r+dr, you'd expect to find about rho * V particles inside, where rho 
          is the number density and V is the volume of that shell. 
  """
  
  distance_data = get_all_distances_from_to(trajectory, from_element, to_element)
  
  interval_width = 0.05
  maximum_distance = sorted(trajectory.cell.cell_dimensions)[0]/2.0 # Smallest dimension #anothercellassumption
  minimum_distance = 1.5  #This'll do for now - it's in Bohr, so that's pretty small
  
  intervals = range(0, int((maximum_distance-minimum_distance)/interval_width))
  interval_rs = list((minimum_distance+(x * interval_width) for x in intervals))
  number_of_intervals = len(intervals)
  rdf = [0] * number_of_intervals # We're not going to use this yet, but we
                                  #  need a list containing the appropriate
                                  #  number of zeros
  
  f = lambda x,minx,maxx: minx<=x<maxx
  
  #particle_coeff = 1 / (
  #                   len(filter(lambda x:x.element==from_element, trajectory.atoms.values() )) * #step 2
  #                   len(trajectory.atoms) / (trajectory.cell.a*trajectory.cell.b*trajectory.cell.c) #step 4
  #                   )
  
  particle_coeff = 1 / (
                     len(filter(lambda x:x.element==from_element, trajectory.atoms.values() )) ** 2 #step 2
                     / (trajectory.cell.a*trajectory.cell.b*trajectory.cell.c) #step 4
                     )
  
  rdf=dict()
  r=interval_width
  while r<maximum_distance:
    atoms_in_interval = 0
    for set in distance_data:
      for distance in set.distances:
        if r<=distance<(r+interval_width):
          atoms_in_interval = atoms_in_interval+1
    rdf[r] = atoms_in_interval * particle_coeff / (trajectory.number_of_steps()*4*numpy.pi*(interval_width**3))
    r=round(r+interval_width,5)
  
  for key in sorted(rdf.keys()):
    print("%f %f" % (key, rdf[key]))
  
    
  

def create_trajectory_from_md_file(filename):
  """
  Creates a trajectory object from a provided md filename.
  
  The format goes, for example for water 
  (spaces replaced by underscores for clarity):
_BEGIN_header
__
_END_header
__
______________________0.0000000000000000E+000
_____________________-1.7384544767983016E+001___-1.7381694634353554E+001____2.8501336294601484E-003__<--_E
______________________9.5004454315338277E-004________________________________________________________<--_T
______________________3.0235618153336794E+001____0.0000000000000000E+000____0.0000000000000000E+000__<--_h
______________________0.0000000000000000E+000____3.0235618153336794E+001____0.0000000000000000E+000__<--_h
______________________0.0000000000000000E+000____0.0000000000000000E+000____3.0235618153336794E+001__<--_h
_H_______________1____2.3621576682294371E-001___-2.1051549139260746E+000____1.1300562284809627E-001__<--_R
_H_______________2___-2.2002081384956260E+000___-3.2233058677591608E+000____1.2188733568063896E+000__<--_R
_O_______________1___-1.5932281040673903E+000___-1.9910154553972281E+000____1.8897261345835500E-004__<--_R
_H_______________1___-8.2511691729044638E-004___-3.6013169146035726E-004___-9.7523996017269724E-004__<--_V
_H_______________2____5.9150969792237012E-004____9.3982676102281318E-004____2.1622400926196277E-004__<--_V
_O_______________1____1.4716930678016571E-005___-3.6519985025362337E-005____4.7816951733250353E-005__<--_V
_H_______________1____4.0287014812819683E-004____7.3424151330816302E-006___-7.3424151330816302E-006__<--_F
_H_______________2___-2.2042797182425522E-004___-2.7658047903899658E-004____2.7658047903899658E-004__<--_F
_O_______________1___-1.8244217630394155E-004____2.6923806390591494E-004___-2.6923806390591494E-004__<--_F
__

The first line is the current time, E is TotalE,HamiltE,KE
T is temperature, h are the cell vectors, R positions, V velocities, 
and F forces.

Records are separated by a blank line with two spaces on, and the final line 
is the same blank line with two spaces.

Atoms are sorted.
  """
  
  trajectory = Trajectory()
  
  try:
    file = open(filename, "r")
  except:
    sys.stderr.write("Could not open file \""+filename+"\" for reading.\n")
    sys.exit(5)
  
  try:
    lines = file.readlines()
  except:
    sys.stderr.write("Could not read file into memory. Maybe it's too big.\n")
    sys.exit(5)

  file.close()


  line_iterator = lines.__iter__()


  #There isn't always a header, apparently.
  line = line_iterator.next()
  if line == " BEGIN header\n":
    # Skip over the header.
    while line != " END header\n":
      line = line_iterator.next()
    line = line_iterator.next() # (And skip the blank line after)
    line = line_iterator.next() # (And this gets us the first data line)
  
  # We need to create some things on the first step, so we do it separately.
    
  trajectory.timesteps.append( float(line.split()[0]))
  
  (h,i,j,arrow,type) = line_iterator.next().split()
  trajectory.energies.append( [ float(h), float(i), float(j) ] )
  
  trajectory.temperatures.append(float(line_iterator.next().split()[0]))
  
  # I am making another cell assumption here, namely that the cell vectors are
  #   coincident with the x y and z vectors. #anothercellassumption
  if trajectory.cell == None:
    trajectory.cell = Cell(float(line_iterator.next().split()[0]),
                           float(line_iterator.next().split()[1]),
                           float(line_iterator.next().split()[2]), 90, 90, 90)
  
  # Now make the atoms and give them their first position.
  (element,index,x,y,z,arrow,type) = line_iterator.next().split()
  while type == "R":
    label = element + " " + index
    trajectory.atoms[label] = MD_Atom(element, index)
    
    trajectory.atoms[label].positions.append(coord_mk(x,y,z))                                           
    (element,index,x,y,z,arrow,type) = line_iterator.next().split()
  
  while type == "V":
    label = element + " " + index
    trajectory.atoms[label].velocities.append(coord_mk(x,y,z))
    (element,index,x,y,z,arrow,type) = line_iterator.next().split()
    
  while type == "F":
    label = element + " " + index
    trajectory.atoms[label].velocities.append(coord_mk(x,y,z))
    
    line = line_iterator.next()
    if line != "  \n":
      (element,index,x,y,z,arrow,type) = line.split()
    else: 
      break
  
  # Now we've done the first step, we can just repeat until we run out.
  for line in line_iterator:
    trajectory.timesteps.append( float(line.split()[0]))
    
    (h,i,j,arrow,type) = line_iterator.next().split()
    trajectory.energies.append( [ float(h), float(i), float(j) ] )
    
    trajectory.temperatures.append(float(line_iterator.next().split()[0]))
  
    # Skip the cell def
    line_iterator.next()
    line_iterator.next()
    line_iterator.next()
    
    (element,index,x,y,z,arrow,type) = line_iterator.next().split()
    while type == "R":
      label = element + " " + index
      trajectory.atoms[label].positions.append(coord_mk(x,y,z))
      (element,index,x,y,z,arrow,type) = line_iterator.next().split()
      
    while type == "V":
      label = element + " " + index
      trajectory.atoms[label].velocities.append(coord_mk(x,y,z))
      (element,index,x,y,z,arrow,type) = line_iterator.next().split()
    
    while type == "F":
      label = element + " " + index
      trajectory.atoms[label].velocities.append(coord_mk(x,y,z))
    
      line = line_iterator.next()
      if line != "  \n":
        (element,index,x,y,z,arrow,type) = line.split()
      else:
        break
    
  # And that should be the end of the file.
  #print(repr(trajectory.atoms))
  return trajectory

if __name__=="__main__":
  sys.stderr.write("Creating trajectory object...\n")
  trajectory = create_trajectory_from_md_file(sys.argv[1])
  sys.stderr.write("Getting distance data...\n")
  distance_data = get_all_distances_from_to(trajectory, "Pu", "O")
  sys.stderr.write("Making graph data files...\n")
  graph_distances(distance_data)
  graph_energies(trajectory)
  graph_temperatures(trajectory)
  #calculate_rdf(trajectory, "Pu", "O")
  sys.stderr.write("Done.\n")
  
  
  


