# thamesmix 0.1.2

This patch is set quite early (less than 1 week after the package has been 
accepted). Unfortunately, the errors found are quite important, which is why I 
sent the patch so soon.

## Patch
This is a patch. In this version I have removed the following bugs:

  -   the center of the ellipsoid was not calculated correctly in some cases, 
      leading to it being empty (now it is calculated correctly)
  -   the permutation problem took longer to solve than it should 
      (now it is computed more quickly)
  -   the adjusted radius was not transferred to the thamesmix function, 
      making the ellipse too large in some cases (now the proper radius is used)
