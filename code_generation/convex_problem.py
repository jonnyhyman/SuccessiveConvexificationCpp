import cvx_sym as cvx

class Definition:

    def __init__(self):

        K = 50

        self.Xv       = cvx.Variable((K, 14), name = 'Xv')
        self.Uv       = cvx.Variable((K, 14), name = 'Uv' )
        self.nuv      = cvx.Variable((14 * (K - 1)), name = 'nuv')
        self.delta    = cvx.Variable(K, name = 'delta')
        self.sigmav   = cvx.Variable(name = 'sigmav')
        self.delta_sv = cvx.Variable(name = 'delta_sv')

        self.x_init   = cvx.Parameter((1, 14), name = 'x_init')
        self.x_final  = cvx.Parameter((1, 14), name = 'x_final')

        # Indexes to expect config variables on
        m_dry = 0
        tan_gamma_gs = 1
        cos_theta_max = 2
        cos_delta_max = 3
        w_B_max = 4
        T_min = 5
        T_max = 6

        self.config   = cvx.Parameter((7), name = 'config')

        # Parameters:
        self.A_bar_parm = cvx.Parameter((K, 14 * 14), name = 'A_bar_parm')
        self.B_bar_parm = cvx.Parameter((K, 14 * 3), name = 'B_bar_parm')
        self.C_bar_parm = cvx.Parameter((K, 14 * 3), name = 'C_bar_parm')
        self.S_bar_parm = cvx.Parameter((K, 14), name = 'S_bar_parm')
        self.z_bar_parm = cvx.Parameter((K, 14), name = 'z_bar_parm')

        self.X_last_parm        = cvx.Parameter((K, 14), name = 'X_last_parm')
        self.U_last_parm        = cvx.Parameter((K, 14), name = 'U_last_parm')
        self.s_last_parm        = cvx.Parameter(name = 's_last_parm')
        self.w_delta_parm       = cvx.Parameter(name = 'w_delta_parm')
        self.w_nu_parm          = cvx.Parameter(name = 'w_nu_parm')
        self.w_delta_sigma_parm = cvx.Parameter(name = 'w_delta_sigma_parm')

        """ The start/end indexes of each variable in the vector V = XABCSz """
        self.idx  = [ 14 ]                # end of x (14,1)
        self.idx += [self.idx[0] + (14 * 14)]  # end of A (14,14)
        self.idx += [self.idx[1] + (14 * 3)]   # end of B (14,3)
        self.idx += [self.idx[2] + (14 * 3)]   # end of C (14,3)
        self.idx += [self.idx[3] + (14 * 1)]   # end of S (14,1)
        self.idx += [self.idx[4] + (14 * 1)]   # end of z (14,1)

        self.constraints = []

        # ----------------------------------------------------- Dynamics:
        for k in range(K - 1):
            self.constraints += [
                self.Xv[k + 1, :] == (
                cvx.reshape(self.A_bar_parm[k, :], (14,14)) * self.Xv[k, :]
                + cvx.reshape(self.B_bar_parm[k, :], (14, 3)) * self.Uv[k, 0:3]
                + cvx.reshape(self.C_bar_parm[k, :], (14, 3)) * self.Uv[k + 1, 0:3]
                + self.S_bar_parm[k, :] * self.sigmav
                + self.z_bar_parm[k, :]
                + self.nuv[k*14 : (k + 1)*14]
                )
            ]

        # ----------------------------------------------------- State constraints:
        self.constraints += [self.Xv[:, 0] >= self.config[m_dry]]

        for k in range(K):
            self.constraints += [
                self.config[tan_gamma_gs] * cvx.norm(self.Xv[k, 2: 4]) <= self.Xv[k, 1],
                self.config[cos_theta_max] <= 1 - 2 * cvx.sum_squares(self.Xv[k, 9:11]),
                cvx.norm(self.Xv[k, 11: 14]) <= self.config[w_B_max]
            ]

            # ------------------------------------------------- Control constraints:

            #B_g = self.U_last_parm[k, 0:3] / cvx.norm(self.U_last_parm[k, 0:3])

            self.constraints += [
            #    self.config[T_min] <= B_g * self.Uv[k, 0:3].T,
                self.Uv[k,0] >= 0,
                cvx.norm(self.Uv[k, 0:3]) <= self.config[T_max],
                self.config[cos_delta_max] * cvx.norm(self.Uv[k, 0:3]) <= self.Uv[k, 0]
            ]

            # ------------------------------------------------- Trust regions:

            dx = self.Xv[k, :] - self.X_last_parm[k, :]
            du = self.Uv[k, :] - self.U_last_parm[k, :]

            self.constraints += [
                cvx.sum_squares(dx) + cvx.sum_squares(du) <= self.delta[k]
            ]

        self.constraints += [cvx.norm(self.sigmav - self.s_last_parm) <= self.delta_sv]
        self.constraints += [self.sigmav >= 0]

        self.constraints += [
            self.Xv[0, 0:7]    == self.x_init[0,0:7],   # init all but attitude
            self.Xv[0, 11:14]  == self.x_init[0,11:14],

            self.Xv[K-1, 1:14] == self.x_final[0,1:14],   # final all but mass
            self.Uv[K-1, 1] == 0,
            self.Uv[K-1, 2] == 0,
        ]

        objective = cvx.Minimize(
            self.sigmav  # minimize flight time
            #-Xv[-1,0]  # minimize fuel use

                # virtual control
                + self.w_nu_parm     * cvx.norm(self.nuv, kind = 1)

                 # trust region on dynamics
                + self.w_delta_parm  * cvx.norm(self.delta)

                # trust region on sigma
                + self.w_delta_sigma_parm * cvx.norm(self.delta_sv, kind = 1)
        )

        self.problem = cvx.Problem(objective, self.constraints)

        gen = cvx.Generate(
                self.problem,
                name = 'Aeroland_core',
                folder = 'raw',
                verbose = 2  # show each stage of the process
              )

Definition()
