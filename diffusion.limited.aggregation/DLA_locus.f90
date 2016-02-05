! Program to simulate a Diffusion limited aggregation (DLA)
! source: http://paulbourke.net/fractals/dla/
!
! Input file format:
!-------
! No of particles
! Size of matrix
! No of stickiness parameters
! Stickiniess parameter 1
! Stickiniess parameter 2
! ...
!
! Output:
!--------
! !) Matrix with integers 1 or 0. Output is ASCII text file
!    File name - fort.3?
! 2) Vector with estimate of stickiness parameter along with vectors used for calculating the estimate.
!    File name - fort.10
!
! gnuplot (http://www.gnuplot.info/) will be used to visualize the binary output file.
! > plot 'fort.3?' matrix with image
!**************************************************************************************************

program dla
implicit none
integer :: input_rows, n_rows, n_cols, mid_row, mid_col, row, col, row_temp, col_temp
integer :: side_of_matrix, rn_1_to_4, n_particles
integer :: temp, no_of_stick_param, idx_stick, idx_particles
integer :: idx_1, idx_2, count_particles
integer :: ii, n, clock
integer, allocatable :: seed(:)
integer(kind=8) :: irec
real, allocatable :: found_neibr_count(:), accepted_neibr_count(:), stick_estimate(:)
real, allocatable :: stickiness(:)
integer(kind=2), allocatable :: matrix(:, :)
real(kind=8):: i
logical :: got_empty_cell, side_1_full, side_2_full, side_3_full, side_4_full
logical :: consider_moving, check(4) ! 4 directions to check
character(len=10) :: fmt
character(len=100)::file1


! seeding the random number generator
call random_seed(size = n)
allocate(seed(n))
call system_clock(count=clock)
seed = clock + 37 * (/ (ii - 1, ii = 1, n) /)
call random_seed(put = seed)

! opening and reading from the input file
call getarg(1,file1)
open(unit=1,file=file1,action='read') 
read(1,*)n_particles
read(1,*)input_rows
read(1,*)no_of_stick_param
allocate(found_neibr_count(no_of_stick_param))
allocate(accepted_neibr_count(no_of_stick_param))
allocate(stick_estimate(no_of_stick_param))
allocate(stickiness(no_of_stick_param))
do ii = 1, no_of_stick_param
    read(1,*)stickiness(ii)
enddo
if (mod(input_rows, 2) == 0) then
    n_rows = input_rows + 1        ! using odd number of rows and cols
else
    n_rows = input_rows
endif
n_cols = n_rows 
mid_row = (n_rows / 2) + 1
mid_col = mid_row
allocate(matrix(n_rows, n_cols))
write (fmt,'(a,i6,a)') '(',n_cols,'i2)'

found_neibr_count = 0.0  ! no. of neighbouring 1s encountered (no of stickable cells)
accepted_neibr_count = 0.0 ! no. of accepted stcking points
do idx_stick = 1, no_of_stick_param
    matrix = 0 ! initialize matrix to 0
    matrix(mid_row, mid_col) = 1 ! center of matrix is 1
    count_particles = 0
    do idx_particles = 1, n_particles
        if (mod(idx_particles, 1000) == 0) then
            write(*,*) idx_stick, idx_particles  ! out the progress of execution of code
        endif
        consider_moving = .true.
        do ! choose side
            !-------------------------------------------------
            ! Choose a side of matrix and add a position
            ! Particle is created on south-side of matrix when random numb. is 1 
            ! 2 - west, 3 - north and 4 - east
            got_empty_cell = .false.
            side_1_full = .false.
            side_2_full = .false.
            side_3_full = .false.
            side_4_full = .false.
            call random_number(harvest = i)
            rn_1_to_4 = int(i * 4) + 1  ! random number from {1, 2, 3 and 4}
            side_of_matrix = rn_1_to_4
            if (side_of_matrix == 1) then
                do  ! choose a position and check if it is already occupied
                    row_temp = n_rows
                    col_temp = rn_1_to_grid(n_cols) 
                    if (matrix(row_temp, col_temp) == 0) then
                        row = row_temp
                        col = col_temp
                        matrix(row, col) = 1
                        got_empty_cell = .true.
                        exit
                    elseif (all(matrix(row_temp, :) == 1)) then
                        side_1_full = .true.
                        exit
                    endif
                enddo
            elseif (side_of_matrix == 2) then
                do
                    row_temp = rn_1_to_grid(n_rows) 
                    col_temp = 1
                    if (matrix(row_temp, col_temp) == 0) then
                        row = row_temp
                        col = col_temp
                        matrix(row, col) = 1
                        got_empty_cell = .true.
                        exit
                    elseif (all(matrix(:, col_temp) == 1)) then
                        side_2_full = .true.
                        exit
                    endif
                enddo
            elseif (side_of_matrix == 3) then
                do 
                    row_temp = 1
                    col_temp = rn_1_to_grid(n_cols) 
                    if (matrix(row_temp, col_temp) == 0) then
                        row = row_temp
                        col = col_temp
                        matrix(row, col) = 1
                        got_empty_cell = .true.
                        exit
                    elseif (all(matrix(row_temp, :) == 1)) then
                        side_3_full = .true.
                        exit
                    endif
                enddo
            else
                do
                    row_temp = rn_1_to_grid(n_rows)
                    col_temp = n_cols
                    if (matrix(row_temp, col_temp) == 0) then
                        row = row_temp
                        col = col_temp
                        matrix(row, col) = 1
                        got_empty_cell = .true.
                        exit
                    elseif (all(matrix(:, col_temp) == 1)) then
                        side_4_full = .true.
                        exit
                    endif
                enddo
            endif
            if (got_empty_cell .eqv. .true.) then
                exit
            elseif (side_1_full .or. side_2_full .or. side_3_full .or. side_4_full .eqv. .true.) then
                consider_moving = .false.
                exit
            endif
        enddo
        !------------------------------------------------------------
        ! 2D random walk
        ! Particles bounce back randomly at the boundaries. No periodic boundary condition.
        if (consider_moving .eqv. .true.) then
            count_particles = count_particles + 1
            do 
                call is_neighbour_filled(matrix, row, col, n_rows, n_cols, check)
                if (any(check(:) .eqv. .true.)) then 
                    found_neibr_count(idx_stick) = found_neibr_count(idx_stick) + 1
                    call random_number(harvest = i)
                    if (i < stickiness(idx_stick)) then
                        accepted_neibr_count(idx_stick) = accepted_neibr_count(idx_stick) + 1
                        exit
                    else
                        call move_particle(matrix, row, col, check, n_rows, n_cols)
                    endif
                else 
                    call move_particle(matrix, row, col, check, n_rows, n_cols)
                endif
            enddo
        endif
    enddo
    do idx_1 = 1, n_rows
        write(30+idx_stick, fmt) (matrix(idx_1, idx_2), idx_2 = 1, n_cols)
    enddo
enddo

stick_estimate(:) = accepted_neibr_count(:) / found_neibr_count(:)
write(10,*)accepted_neibr_count
write(10,*)found_neibr_count
write(10,*)stick_estimate
write(10,*)count_particles

contains


!------------------------------------------------------------------
! subroutine to check if neighbour cell is filled with 1
!---------------------------------------------------------------------
subroutine is_neighbour_filled(matrix, row, col, n_rows, n_cols, check)
implicit none
integer(kind=2), intent(in), allocatable :: matrix(:, :)
integer, intent(in) :: row, col, n_rows, n_cols
logical, intent(out) :: check(4) ! 4 directions to check in square lattice
integer :: i, row_temp, col_temp

check = .false. ! default is "neighbour cell is NOT filled"
do i = 1, 4 ! 4 is the no. of directions possible in a square lattice
    if(i == 1) then ! south
        row_temp = row + 1
        if (row_temp <= n_rows) then
            if(matrix(row_temp, col) == 1) then
                check(i) = .true.
            endif
        endif
    endif
    if(i == 2) then ! west
        col_temp = col - 1
        if (col_temp >= 1) then
            if(matrix(row, col_temp) == 1) then
                check(i) = .true.
            endif
        endif
    endif
    if(i == 3) then ! north
        row_temp = row - 1
        if (row_temp >= 1) then
            if(matrix(row_temp, col) == 1) then
                check(i) = .true.
            endif
        endif
    endif
    if(i == 4) then ! east
        col_temp = col + 1
        if (col_temp <= n_cols) then
            if(matrix(row, col_temp) == 1) then
                check(i) = .true.
            endif
        endif
    endif
enddo
end subroutine is_neighbour_filled


!---------------------------------------------------------------
! subroutine to move the particle
! direction: 1 - south, 2 - west, 3 - north, 4 - east. 
!---------------------------------------------------------------
subroutine move_particle(matrix, row, col, check, n_rows, n_cols)
implicit none
integer(kind=2), intent(inout), allocatable :: matrix(:, :)
integer, intent(inout) :: row, col
integer, intent(in) :: n_rows, n_cols
logical, intent(in) :: check(4) ! true if neighbouring cell is filled with 1
real(kind=8) :: i
integer :: direction, row_temp, col_temp

do
    call random_number(harvest = i)
    direction = int(i * 4) + 1  
    if (check(direction) .eqv. .false.) then
        exit
    endif
enddo

if (direction == 1) then ! move south
    if (row < n_rows) then  
        row_temp = row + 1
        matrix(row_temp, col) = 1  ! new cell now has 1 
        matrix(row, col) = 0       ! old cell now has 0
        row = row_temp             ! update cell
    elseif (row == n_rows) then 
       ! optional boundary condition
    else
        write(*,*)"error: row > n_rows"
        STOP
    endif
elseif (direction == 2) then
    if (col > 1) then
        col_temp = col - 1
        matrix(row, col_temp) = 1
        matrix(row, col) = 0
        col = col_temp
    elseif (col == 1) then
       ! optional boundary condition
    else
        write(*,*)"error: col < 1"
        STOP
    endif
elseif (direction == 3) then
    if (row > 1) then  
        row_temp = row - 1
        matrix(row_temp, col) = 1
        matrix(row, col) = 0
        row = row_temp
    elseif (row == 1) then
       ! optional boundary condition
    else
        write(*,*)"error: row < 1"
        STOP
    endif
else
    if (col < n_cols) then
        col_temp = col + 1
        matrix(row, col_temp) = 1
        matrix(row, col) = 0
        col = col_temp
    elseif (col == n_cols) then
        ! optional boundary condition
    else
        write(*,*)"error: col > n_cols"
        STOP
    endif
endif

end subroutine move_particle


!---------------------------------------
! function to generate random numbers between 1 and grid_length included
!------------------------------------------
integer function rn_1_to_grid(grid_length)
implicit none
integer, intent(in) :: grid_length
real(kind=8):: i
call random_number(harvest = i)
rn_1_to_grid = int(i * grid_length) + 1 
end function rn_1_to_grid 

end program dla
