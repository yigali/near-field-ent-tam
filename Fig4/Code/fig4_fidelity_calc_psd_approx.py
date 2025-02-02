


# packages
import torch

# cuda init
cuda0 = torch.device('cuda:0')


def get_matrix_square_root(input_mat):

    # this works only for PSD matrices
    # the linalg.eigh is more stable than the linald.svd since svd becomes unstable for small eigenvalues
    # both however become unstable when there are many repeated eigenvalues
    
    regularization_matrix = 1e-10 * torch.diag(torch.rand(input_mat.size(0), dtype=torch.float64, device=cuda0))
    
    # get the decomposition
    L, Q = torch.linalg.eigh(input_mat + regularization_matrix)

    # get the square root of eigenvalues, discarding any numerically zero but negative values
    # this holds as long as matrices in the loss function are valid density matrices
    L_sqrt = torch.where(L > 0, torch.sqrt(L.to(dtype=torch.float64)), 0.)
        
    # use the decomposition to build back the matrix
    output_mat = torch.matmul(Q, torch.matmul(torch.diag(L_sqrt), Q.T.conj()))
    
    
    # # get the decomposition
    # U, S, V = torch.linalg.svd(input_mat + regularization_matrix)

    # # get the square root of eigenvalues, discarding any numerically zero but negative values
    # # this holds as long as matrices in the loss function are valid density matrices
    # L_sqrt = torch.where(S > 0, torch.sqrt(S.to(dtype=torch.complex128)), 0.)
        
    # # use the decomposition to build back the matrix
    # output_mat = torch.matmul(U, torch.matmul(torch.diag(L_sqrt), V))


    return output_mat

def fidelity_loss(rho_true_mat, rho_est_mat):
    
    # this is loss based on fidelity between true and estimated density matrices (note this is not symmetric)
    
    # rho_est_mat = op_to_mat_torch(rho_est)
    # rho_true_mat = op_to_mat_torch(rho_true)
    
    # get the square root of rho_true
    sqrt_rho_true_mat = get_matrix_square_root(rho_true_mat)
    
    # take three matrix products sqrt_rho_true * rho_est * sqrt_rho_true
    product_1 = torch.matmul(sqrt_rho_true_mat, rho_est_mat)
    product_2 = torch.matmul(product_1, sqrt_rho_true_mat)


    # take sqrtm of the product
    sqrtm_fidelity_arg = get_matrix_square_root(product_2)
    
    # take trace and square it
    fidelity_val = torch.trace(sqrtm_fidelity_arg)
    
    return fidelity_val.real ** 2



def is_psd_and_hermitian(matrix):
    # Check if the matrix is real and square
    if matrix.shape != (2, 2) or not torch.is_tensor(matrix):
        raise ValueError("Input must be a 2x2 tensor")

    # Check if the matrix is Hermitian
    if not torch.equal(matrix, matrix.T):
        return False, "Matrix is not Hermitian"

    # Calculate eigenvalues
    eigenvalues, _ = torch.linalg.eigh(matrix)

    # Check if all eigenvalues are non-negative
    if torch.all(eigenvalues >= -1e-4):
        return True, "Matrix is PSD and Hermitian"
    else:
        return False, "Matrix is Hermitian but not PSD"




if __name__ == '__main__':

    
    # ideal mat
    ideal_mat = 0.5 * torch.tensor([[1, -1], [-1, 1]], dtype=torch.float64, device=cuda0)
    
    # tomography mat
    tomography_mat = torch.tensor([[0.9553, -0.9292], [-0.9428, 1]], dtype=torch.float64, device=cuda0)
    tomography_mat[1, 0] = (1/2) * (-0.9292 + -0.9428)
    tomography_mat[0, 1] = tomography_mat[1, 0]
    
    
    # example num of bootstrap iterations
    bootstrap_num = 1000
    
    
    # initialize purity vector and fidelity vector
    fidelity_vector = torch.zeros(bootstrap_num, dtype=torch.float64, device=cuda0)
    purity_vector = torch.zeros(bootstrap_num, dtype=torch.float64, device=cuda0)
    
    
    # start the bootstrap loop
    for bootstrap_step in range(bootstrap_num):

        # joint error
        joint_std = (1/2) * (0.0402 + 0.0859)
    
        # error mat elements
        # error_mat[0, 0] = torch.normal(torch.zeros(1, dtype=torch.float64, device=cuda0), 0.0103 * torch.ones(1, dtype=torch.float64, device=cuda0))
        
        error_mat = torch.zeros(tomography_mat.size(), dtype=torch.float64, device=cuda0)
        error_mat[0, 0] = -0.0103 + 2 * 0.0103 * torch.rand(1, dtype=torch.float64, device=cuda0)
        error_mat[0, 1] = -joint_std + 2 * joint_std * torch.rand(1, dtype=torch.float64, device=cuda0)
        error_mat[1, 0] = error_mat[0, 1]
        error_mat[1, 1] = -0.0224 + 2 * 0.0224 * torch.rand(1, dtype=torch.float64, device=cuda0)
        
        
        # get the bootstrapped tomography mat
        tomography_mat = tomography_mat + error_mat
        
        # normalize to trace 1
        tomography_mat = tomography_mat / torch.trace(tomography_mat)
        
        
        # check if the matrix is PSD and Hermitian
        psd_hermitian_test, test_message = is_psd_and_hermitian(tomography_mat)
        print(test_message, "at bootstrap step", bootstrap_step)
        
        # if not PSD then do the optimization to find
        # the closest PSD approximation to the current realization
        if psd_hermitian_test == False:
            
            # do several trials for noiseless
            current_best_final = 100
            current_best_mat_final = torch.zeros(tomography_mat.size(), dtype=torch.complex128, device=cuda0)
            num_of_trials = 1
            
            for trial_id in range(0, num_of_trials, 1):
            

                # variational mat
                # estimated_mat = (-1 + 2 * torch.rand(tomography_mat.size(), dtype=torch.float64, device=cuda0)) +\
                    # 1j * (-1 + 2 * torch.rand(tomography_mat.size(), dtype=torch.float64, device=cuda0))
                    
                estimated_mat = -1 + 2 * torch.rand(tomography_mat.size(), dtype=torch.float64, device=cuda0)
                    
            
            
                # turn on the gradient
                estimated_mat.requires_grad = True
            
                # Define optimizer and learning rates
                scheduler_steps=[1500, 2500]
                scheduler_gamma = 0.15
                learning_rate = 1e-2
                optimizer = torch.optim.Adam([estimated_mat], lr = learning_rate)
                scheduler = torch.optim.lr_scheduler.StepLR(optimizer, 1, gamma=scheduler_gamma)
            
                
                # best coeffs and best loss save
                best_loss = 100
                best_step_id = 0
                best_estimated_mat = torch.zeros(estimated_mat.size(), dtype=torch.complex128, device=cuda0)
            
                # example ptimization steps number
                optim_steps = 3000
                
                
                # optimizer loop
                for step_id in range(optim_steps):
                    
                    # zero the gradients
                    optimizer.zero_grad()
                    
                    # make the estimated matrix PSD and Hermitian and trace one
                    estimated_mat_regularized = torch.matmul(torch.conj(torch.permute(estimated_mat, dims = (1, 0))), estimated_mat).real
                    estimated_mat_regularized = estimated_mat_regularized / torch.trace(estimated_mat_regularized)
            
                    # get the loss as Frobenius norm
                    loss = 0.5 * torch.linalg.matrix_norm(estimated_mat_regularized - tomography_mat, "nuc")
                    # print(torch.linalg.matrix_norm(estimated_mat_regularized.detach() - ideal_mat))
                    
            
                    # save the lowest loss argument 
                    if loss.tolist() < best_loss:
                        
                        best_loss_step = step_id
                        best_loss = loss.tolist()
                        best_estimated_mat = estimated_mat.detach()
                        
                    
                    # get grads and step
                    loss.backward()
                    optimizer.step()
                    
                    
                    if step_id % 1000 == 0:
                    
                        print("Loss is", loss.tolist(), "at step", step_id, "trial", trial_id)
            
                    
                       
                    # print stats about the learning rate update
                    if scheduler_gamma is not None and step_id in scheduler_steps:
                        scheduler.step()
                        print(f"lr={optimizer.state_dict()['param_groups'][0]['lr']}")
                        

                # check best so far
                if current_best_final > best_loss:
                    
                    current_best_mat_final = best_estimated_mat
                     
    
            # now take the best guess, make it PSD and Hermitian and trace 1 and get the fidelity to the ideal mat
            print(current_best_mat_final)
            estimated_mat_final = torch.matmul(torch.conj(torch.permute(current_best_mat_final, dims = (1, 0))), current_best_mat_final).real
            estimated_mat_final = estimated_mat_final / torch.trace(estimated_mat_final)
            
            # get the fidelity and purity between the best guess and the ideal mat
            current_fidelity = fidelity_loss(ideal_mat, estimated_mat_final)
            current_purity = torch.trace(torch.matmul(estimated_mat_final, estimated_mat_final)).real
            
            # save them into the tensor
            fidelity_vector[bootstrap_step] = current_fidelity
            purity_vector[bootstrap_step] = current_purity
            
            # save on disk
            torch.save(fidelity_vector, 'fidelity_vector.pt')
            torch.save(purity_vector, 'purity_vector.pt')
            
            # print stats for the current bootstrap step
            print("\nFidelity", current_fidelity.tolist(), "at bootstrap step", bootstrap_step)
            print("Purity", current_purity.tolist(), "at bootstrap step", bootstrap_step, "\n")
            
            
        # if PSD then just get the fidelity and purity
        else:
            
            
            # get the fidelity and purity between the best guess and the ideal mat
            current_fidelity = fidelity_loss(ideal_mat, tomography_mat)
            current_purity = torch.trace(torch.matmul(tomography_mat, tomography_mat)).real
            
            # save them into the tensor
            fidelity_vector[bootstrap_step] = current_fidelity
            purity_vector[bootstrap_step] = current_purity
            
            # save on disk
            torch.save(fidelity_vector, 'fidelity_vector.pt')
            torch.save(purity_vector, 'purity_vector.pt')
            
            # print stats for the current bootstrap step
            print("\nFidelity", current_fidelity.tolist(), "at bootstrap step", bootstrap_step)
            print("Purity", current_purity.tolist(), "at bootstrap step", bootstrap_step, "\n")
        
        























