function CheckUnstructuredMesh2d(node,elem,edge,face)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the grid data using the following tests:
%
% 1. Directed area must sum up to zero around every node.
% 2. Directed area must sum up to zero over the entire grid.
% 3. Global sum of the boundary normal vectors must vanish.
% 4. Global sum of the boundary face normal vectors must vanish.
% 5. Check element volumes which must be positive.
% 6. Check dual volumes which must be positive.
% 7. Global sum of the dual volumes must be equal to the sum of element
%     volumes. 
%
% Add more tests you can think of.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 fprintf('Checking grid data....\n');

  mag_dav = zero
  mag_bn  = zero

%--------------------------------------------------------------------------------
% Directed area sum check
%--------------------------------------------------------------------------------

% Compute the sum of the directed area for each node.

   sum_dav_i = zero
  do i = 1, nedges
   n1 = edge(i)%n1
   n2 = edge(i)%n2
   sum_dav_i(n1,:) = sum_dav_i(n1,:) + edge(i)%dav(:)*edge(i)%da
   sum_dav_i(n2,:) = sum_dav_i(n2,:) - edge(i)%dav(:)*edge(i)%da
   mag_dav = mag_dav + edge(i)%da
  end do
   mag_dav = mag_dav/real(nedges,p2)

% Add contribution from boundary edges.
  do i = 1, nbound
   do j = 1, bound(i)%nbfaces

     n1 = bound(i)%bnode(j)
     n2 = bound(i)%bnode(j+1)

     sum_dav_i(n1,1) = sum_dav_i(n1,1) + half*bound(i)%bfnx(j)*bound(i)%bfn(j)
     sum_dav_i(n1,2) = sum_dav_i(n1,2) + half*bound(i)%bfny(j)*bound(i)%bfn(j)

     sum_dav_i(n2,1) = sum_dav_i(n2,1) + half*bound(i)%bfnx(j)*bound(i)%bfn(j)
     sum_dav_i(n2,2) = sum_dav_i(n2,2) + half*bound(i)%bfny(j)*bound(i)%bfn(j)

   end do
  end do

% Compute also the sum of the boundary normal vector (at nodes).

  sum_bn = 0
  do i = 1, nbound
   do j = 1, bound(i)%nbnodes
     k = bound(i)%bnode(j)
     if (j > 1 .and. k==bound(i)%bnode(1)) cycle %Skip if the last node is equal to the first node).
    sum_bn(1)      = sum_bn(1)      + bound(i)%bnx(j)*bound(i)%bn(j)
    sum_bn(2)      = sum_bn(2)      + bound(i)%bny(j)*bound(i)%bn(j)
    mag_bn = mag_bn + abs(bound(i)%bn(j))
   end do
    mag_bn = mag_bn/real(bound(i)%nbnodes,p2)
  end do

% Global sum of boundary normal vectors must vanish.

%  if (sum_bn(1) > 1.0e-12_p2*mag_bn .and. sum_bn(2) > 1.0e-12_p2*mag_bn) then
%   write(*,*) "--- Global sum of the boundary normal vector:"
%   write(*,'(a19,es10.3)') "    sum of bn_x = ", sum_bn(1)
%   write(*,'(a19,es10.3)') "    sum of bn_y = ", sum_bn(2)
%   write(*,*) "Error: boundary normal vectors do not sum to zero..."
%   stop
%  endif

% Sum of the directed area vectors must vanish at every node.

  do i = 1, nnodes
   if (abs(sum_dav_i(i,1))>1.0e-12_p2*mag_dav .or. abs(sum_dav_i(i,2))>1.0e-12_p2*mag_dav) then
   write(*,'(a11,i5,a7,2es10.3,a9,2es10.3)') &
    " --- node=", i, " (x,y)=", node(i)%x, node(i)%y, " sum_dav=",sum_dav_i(i,:)
   endif
  end do

   write(*,*) "--- Max sum of directed area vector around a node:"
   write(*,*) "  max(sum_dav_i_x) = ", maxval(sum_dav_i(:,1))
   write(*,*) "  max(sum_dav_i_y) = ", maxval(sum_dav_i(:,2))

  if (maxval(abs(sum_dav_i(:,1)))>1.0e-12_p2*mag_dav .or. &
      maxval(abs(sum_dav_i(:,2)))>1.0e-12_p2*mag_dav) then
   write(*,*) "--- Max sum of directed area vector around a node:"
   write(*,*) "  max(sum_dav_i_x) = ", maxval(sum_dav_i(:,1))
   write(*,*) "  max(sum_dav_i_y) = ", maxval(sum_dav_i(:,2))
   write(*,*) "Error: directed area vectors do not sum to zero..."
   stop
  endif

% Of course, the global sum of the directed area vector sum must vanish.
   sum_dav = zero
  do i = 1, nnodes
   sum_dav = sum_dav + sum_dav_i(i,:)
  end do

   write(*,*) "--- Global sum of the directed area vector:"
   write(*,'(a19,es10.3)') "    sum of dav_x = ", sum_dav(1)
   write(*,'(a19,es10.3)') "    sum of dav_y = ", sum_dav(2)

  if (sum_dav(1) > 1.0e-12_p2*mag_dav .and. sum_dav(2) > 1.0e-12_p2*mag_dav) then
   write(*,*) "Error: directed area vectors do not sum globally to zero..."
   write(*,*) "--- Global sum of the directed area vector:"
   write(*,'(a19,es10.3)') "    sum of dav_x = ", sum_dav(1)
   write(*,'(a19,es10.3)') "    sum of dav_y = ", sum_dav(2)
   stop
  endif

%--------------------------------------------------------------------------------
% Global sum check for boundary face vector
%--------------------------------------------------------------------------------
  sum_bfn = 0
  do i = 1, nbound
   do j = 1, bound(i)%nbfaces
     sum_bfn(1) =  sum_bfn(1) + bound(i)%bfnx(j)*bound(i)%bfn(j)
     sum_bfn(2) =  sum_bfn(2) + bound(i)%bfny(j)*bound(i)%bfn(j)
   end do
  end do

   write(*,*) "--- Global sum of the boundary face vector:"
   write(*,'(a19,es10.3)') "    sum of bfn_x = ", sum_bfn(1)
   write(*,'(a19,es10.3)') "    sum of bfn_y = ", sum_bfn(2)

  if (sum_bfn(1) > 1.0e-12_p2*mag_bn .and. sum_bfn(2) > 1.0e-12_p2*mag_bn) then
   write(*,*) "Error: boundary face normals do not sum globally to zero..."
   write(*,*) "--- Global sum of the boundary face normal vector:"
   write(*,'(a19,es10.3)') "    sum of bfn_x = ", sum_bfn(1)
   write(*,'(a19,es10.3)') "    sum of bfn_y = ", sum_bfn(2)
   stop
  endif

%--------------------------------------------------------------------------------
% Volume check
%--------------------------------------------------------------------------------
% (1)Check the element volume: make sure there are no zero or negative volumes

   vol_min =  1.0e+15
   vol_max = -1.0
   vol_ave =  zero

       ierr = 0
   sum_volc = zero
  do i = 1, nelms

      vol_min = min(vol_min,elm(i)%vol)
      vol_max = max(vol_max,elm(i)%vol)
      vol_ave = vol_ave + elm(i)%vol

   sum_volc = sum_volc + elm(i)%vol

   if (elm(i)%vol < zero) then
     write(*,*) "Negative volc=",elm(i)%vol, " elm=",i, " stop..."
     ierr = ierr + 1
   endif

   if (abs(elm(i)%vol) < 1.0e-14_p2) then
     write(*,*) "Vanishing volc=",elm(i)%vol, " elm=",i, " stop..."
     ierr = ierr + 1
   endif

  end do

   vol_ave = vol_ave / real(nelms)

   write(*,*)
   write(*,'(a30,es25.15)') "    minimum element volume = ", vol_min
   write(*,'(a30,es25.15)') "    maximum element volume = ", vol_max
   write(*,'(a30,es25.15)') "    average element volume = ", vol_ave
   write(*,*)

%--------------------------------------------------------------------------------
% (2)Check the dual volume (volume around a node)

   vol_min =  1.0e+15
   vol_max = -1.0
   vol_ave =  zero

      ierr = 0
   sum_vol = zero
  do i = 1, nnodes

      vol_min = min(vol_min,node(i)%vol)
      vol_max = max(vol_max,node(i)%vol)
      vol_ave = vol_ave + node(i)%vol

   sum_vol = sum_vol + node(i)%vol

   if (node(i)%vol < zero) then
     write(*,*) "Negative vol=",node(i)%vol, " node=",i, " stop..."
     ierr = ierr + 1
   endif

   if (abs(node(i)%vol) < 1.0e-14_p2) then
     write(*,*) "Vanishing vol=",node(i)%vol, " node=",i, " stop..."
     ierr = ierr + 1
   endif

  end do

   vol_ave = vol_ave / real(nnodes)

   write(*,*)
   write(*,'(a30,es25.15)') "    minimum dual volume = ", vol_min
   write(*,'(a30,es25.15)') "    maximum dual volume = ", vol_max
   write(*,'(a30,es25.15)') "    average dual volume = ", vol_ave
   write(*,*)


  if (ierr > 0) stop

  if (abs(sum_vol-sum_volc) > 1.0e-08_p2*sum_vol) then
   write(*,*) "--- Global sum of volume: must be the same"
   write(*,'(a19,es10.3)') "    sum of volc = ", sum_volc
   write(*,'(a19,es10.3)') "    sum of vol  = ", sum_vol
   write(*,'(a22,es10.3)') " sum_vol-sum_volc  = ", sum_vol-sum_volc
   write(*,*) "Error: sum of dual volumes and cell volumes do not match..."
   stop
  endif

  call check_skewness_nc
  call compute_ar

  write(*,*)
  write(*,*) "Grid data look good%"

 end subroutine check_grid_data

%*******************************************************************************
%* Skewness computation for edges.
%*******************************************************************************
 subroutine check_skewness_nc
 use edu2d_my_main_data  , only : nedges, edge
 use edu2d_constants     , only : p2, zero

 implicit none
 integer :: i
 real(p2) :: e_dot_n, e_dot_n_min, e_dot_n_max, alpha

     e_dot_n = zero
 e_dot_n_min = 100000.0_p2
 e_dot_n_max =-100000.0_p2

  edges : do i = 1, nedges

   alpha = edge(i)%ev(1)*edge(i)%dav(1) + edge(i)%ev(2)*edge(i)%dav(2)
   e_dot_n     = e_dot_n + abs(alpha)
   e_dot_n_min = min(e_dot_n_min, abs(alpha))
   e_dot_n_max = max(e_dot_n_max, abs(alpha))

  end do edges

  e_dot_n = e_dot_n / real(nedges,p2)

 write(*,*)
 write(*,*) " ------ Skewness check (NC control volume) ----------"
 write(*,*) "   L1(e_dot_n) = ", e_dot_n
 write(*,*) "  Min(e_dot_n) = ", e_dot_n_min
 write(*,*) "  Max(e_dot_n) = ", e_dot_n_max
 write(*,*) " ----------------------------------------------------"

 end subroutine check_skewness_nc


%*******************************************************************************
%* Control volume aspect ratio
%*******************************************************************************
 subroutine compute_ar
 use edu2d_my_main_data  , only : node, nnodes, elm, nelms
 use edu2d_constants     , only : p2, zero, one, half, two

 implicit none
 integer :: i, n1, n2
 real(p2) :: ar, ar_min, ar_max, nnodes_eff

 integer  :: k
 real(p2) :: side_max, side(4), side_mid, side_min, height

% Initialization

  node1 : do i = 1, nnodes
   node(i)%ar     = zero
  end do node1

% Compute element aspect-ratio: longest_side^2 / vol

  elm1: do i = 1, nelms

   side_max = -one
   
   do k = 1, elm(i)%nvtx

     n1 = elm(i)%vtx(k)
    if (k == elm(i)%nvtx) then
     n2 = elm(i)%vtx(1)
    else
     n2 = elm(i)%vtx(k+1)
    endif

     side(k) = sqrt( (node(n2)%x-node(n1)%x)**2 + (node(n2)%y-node(n1)%y)**2 )
    side_max =  max(side_max, side(k))

   end do

   if (elm(i)%nvtx == 3) then

 % AR for triangle:  Ratio of a half of a square with side_max to volume
    elm(i)%ar = (half*side_max**2) / elm(i)%vol

    if     (side(1) >= side(2) .and. side(1) >= side(3)) then

       side_max = side(1)
      if (side(2) >= side(3)) then
       side_mid = side(2); side_min = side(3)
      else
       side_mid = side(3); side_min = side(2)
      endif

    elseif (side(2) >= side(1) .and. side(2) >= side(3)) then

       side_max = side(2)
      if (side(1) >= side(3)) then
       side_mid = side(1); side_min = side(3)
      else
       side_mid = side(3); side_min = side(1)
      endif

    else

       side_max = side(3)
      if (side(1) >= side(2)) then
       side_mid = side(1); side_min = side(2)
      else
       side_mid = side(2); side_min = side(1)
      endif

    endif

       height = two*elm(i)%vol/side_mid
    elm(i)%ar = side_mid/height

   else

  % AR for quad: Ratio of a square with side_max to volume
    elm(i)%ar = side_max**2 / elm(i)%vol

   endif

  end do elm1

% Compute the aspect ratio:
  node2 : do i = 1, nnodes

    node(i)%ar = zero
   do k = 1, node(i)%nelms
    node(i)%ar = node(i)%ar + elm(node(i)%elm(k))%ar
   end do

    node(i)%ar = node(i)%ar / real(node(i)%nelms, p2)

  end do node2

% Compute the min/max and L1 of AR

  nnodes_eff= zero
         ar = zero
     ar_min = 100000.0_p2
     ar_max =-100000.0_p2

  node3: do i = 1, nnodes
   if (node(i)%bmark /= 0) cycle node3
   ar     = ar + abs(node(i)%ar)
   ar_min = min(ar_min, abs(node(i)%ar))
   ar_max = max(ar_max, abs(node(i)%ar))
   nnodes_eff = nnodes_eff + one
  end do node3

  ar = ar / nnodes_eff

 write(*,*)
 write(*,*) " ------ Aspect ratio check (NC control volume) ----------"
 write(*,*) " Interior nodes only"
 write(*,*) "   L1(AR) = ", ar
 write(*,*) "  Min(AR) = ", ar_min
 write(*,*) "  Max(AR) = ", ar_max

  nnodes_eff= zero
         ar = zero
     ar_min = 100000.0_p2
     ar_max =-100000.0_p2

  node4: do i = 1, nnodes
   if (node(i)%bmark == 0) cycle node4
   ar     = ar + abs(node(i)%ar)
   ar_min = min(ar_min, abs(node(i)%ar))
   ar_max = max(ar_max, abs(node(i)%ar))
   nnodes_eff = nnodes_eff + one
  end do node4

  ar = ar / nnodes_eff

 write(*,*)
 write(*,*) " Boundary nodes only"
 write(*,*) "   L1(AR) = ", ar
 write(*,*) "  Min(AR) = ", ar_min
 write(*,*) "  Max(AR) = ", ar_max
 write(*,*) " --------------------------------------------------------"

 end subroutine compute_ar
