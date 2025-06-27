## R CMD check results

0 errors | 0 warnings | 0 notes

## Patch
This is a patch. In this version I have removed the following bugs:

  -   the center of the ellipsoid was not calculated correctly in some cases, 
      leading to it being empty (now it is calculated correctly)
  -   the permutation problem took longer to solve than it should 
      (now it is computed more quickly)
  -   the adjusted radius was not transferred to the thamesmix function, 
      making the ellipse too large in some cases (now the proper radius is used)
