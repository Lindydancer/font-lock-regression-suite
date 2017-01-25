! Copyright 2005, 2010 Allinea Software Ltd.

! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

program mainf
  include 'mpif.h'

  real :: a(10), b(5, 5), c(13,15,17), conv1(3,3), conv2(3,3)
  integer :: i, j, k, ierr, size, rank, bigdim
  real :: rnum
  real, allocatable :: twod(:,:), threed(:,:,:)
  integer, allocatable :: seed(:)
  real :: t 
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

  if (rank.eq.0) then
    write(*,*) "A real convoluted example code"
  end if

  b(:, :) = 0
  do i = 1, 5
     b(i, i) = 1
  end do
  
  if (rank.eq.0) then
    write (*,*) "Start of input b"
    write (*,*) b
    write (*,*) "End of input b"
  end if

  conv1 = reshape((/ 0, 0, 0, 0, 1, 0, 0, 0, 0 /), shape(conv1))
  conv2 = reshape((/ 0.33, 0.0, 0.0, 0.0, 0.33, 0.0, 0.0, 0.0, 0.33 /), shape(conv2))
  


  
  t = mod(int( mpi_wtime()), 128)

  call random_seed(SIZE = n)
  allocate(seed(n))

  seed = rank + 37 + t;
  
  call random_seed(PUT = seed)
  deallocate(seed)

  call random_number(rnum)

  if (rank.eq.0) then
    write (*,*) rnum;
  end if

  bigdim = rnum * 1200
  allocate(twod(bigdim, bigdim))
  allocate(threed(13,15,17))



  call filltwo(twod, bigdim, bigdim)
  call fillthree(threed, 13, 15, 17)



! Try identity convolution with the identity matrix

  call convolute(b, 5, 5, conv1)

! Did it look good?

  if (rank.eq.0) then
    write(*,*) "Output of the convolution of b"
    write(*,*) b
    write(*,*) "End of the convolution of b"
  end if

! Now lets try the bigger array
  call convolute(twod, bigdim, bigdim, conv1)

  if (rank.eq.0) then
    write(*,*) "Convolution of twod complete"
  end if

  call convolute(twod, bigdim, bigdim, conv2)




  deallocate(twod)
  deallocate(threed)


  if (rank.eq.0) then
    write(*,*) 'This is a test.'
  endif
  call MPI_Finalize()

  stop



end program mainf

subroutine fillone(a, n)
  real :: a(n)
  integer :: i
  do i = 1, n
     a(i) = i
  end do
end subroutine fillone

subroutine filltwo(b, n, m)
  real :: b(n, m)
  integer :: i, j
  do i = 1, n
     do j = 1, m
        b(i, j) = i * m + j
     end do
  end do
end subroutine filltwo

subroutine fillthree(c, n, m, o)
  integer :: n, m, o
  real :: c(n, m, o)
  integer :: i, j, k

  do i = 1, n
     do j = 1, m
        do k = 1, o
           c(i, j, k) = i * m * o + j * o + k
        end do
     end do
  end do
end subroutine fillthree

  

subroutine convolute(b, n, m, conv)
  real :: b(n, m)
  integer :: i, j, k, l
  real :: conv(3, 3)

  real :: c(n, m)

  c(:,:) = 0

  do i = 1, n
     do j = 1, m
        do k = 1, 3
           do l = 1, 3
              c(i, j) = c(i, j) +  &
                 b(i + k - 2, j + l - 2) * conv(k, 1)
           end do
        end do
     end do
  end do

  do i = 1, n
     do j = 1, m
        b(i, j) = c(i, j)
     end do
  end do
end subroutine convolute
