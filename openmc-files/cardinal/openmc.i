pin_pitch = ${fparse 5.5685 / 10 * 1.180}

[Mesh]
  [pin]
    type = PolygonConcentricCircleMeshGenerator
    polygon_size = ${fparse pin_pitch}
    polygon_size_style = radius
    num_sides = 6
    num_sectors_per_side = '2 2 2 2 2 2'
  []
  [bundle]
    type = PatternedHexMeshGenerator
    pattern_boundary = none
    inputs = 'pin'
    pattern = '0 0 0 0 0 0 0 0 0 0;
              0 0 0 0 0 0 0 0 0 0 0;
             0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0;
           0 0 0 0 0 0 0 0 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
           0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0;
             0 0 0 0 0 0 0 0 0 0 0 0;
              0 0 0 0 0 0 0 0 0 0 0;
               0 0 0 0 0 0 0 0 0 0
'
  []
[]

[Problem]
  type = OpenMCCellAverageProblem
  power = 1
  cell_level = 2
  normalize_by_global_tally = false


  [Tallies]
    [cell]
      type = CellTally
      blocks = '1'
      output = 'unrelaxed_tally_rel_error unrelaxed_tally_std_dev'
    []
  []
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
[]
