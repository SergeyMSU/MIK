PROGRAM hello_world_mpi
include 'mpif.h'

integer process_Rank, size_Of_Cluster, ierror, tag
integer :: mas(2, 2), i
integer :: bank(2, 2)
integer status(MPI_STATUS_SIZE)


call MPI_INIT(ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, size_Of_Cluster, ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, process_Rank, ierror)

if(size_Of_Cluster /= 4) STOP "ERROR iuyhk"

mas(1, 1) = 1 + 4 * process_Rank
mas(2, 1) = 2 + 4 * process_Rank
mas(1, 2) = 3 + 4 * process_Rank
mas(2, 2) = 4 + 4 * process_Rank

!print *, 'Hello World from process: ', process_Rank, 'of ', size_Of_Cluster

do i = 1, 3
	if(process_Rank == i) call MPI_SEND(mas, 4, MPI_INTEGER, 0, 100, MPI_COMM_WORLD, ierror)
	if(process_Rank == 0) then
		call MPI_RECV(bank, 4, MPI_INTEGER, i, 100, MPI_COMM_WORLD, status, ierror)
		mas = mas + bank
	end if
end do


if(process_Rank == 0) then
	print*, mas
end if


call MPI_FINALIZE(ierror)
END PROGRAM