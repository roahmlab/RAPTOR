import numpy as np
import math
import pinocchio as pin
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
import scipy.io

def integrate(model, constraint_model, constraint_data,
              ts_sim, x0, 
              ts = None, xs = None, us = None, act_matrix = None,
              Kp = None, Kd = None,
              kind='zero',
              method='RK45'):
    t = 0.0
    nq = model.nq
    nv = model.nv
    data_sim = model.createData()
    pin.initConstraintDynamics(model, data_sim, constraint_model)
    prox_settings = pin.ProximalSettings(1e-12, 1e-12, 100)
    
    def control(t, x):
        u_openloop = interp1d(ts, us.T, kind=kind)(t)
        
        xdes = interp1d(ts, xs.T, kind=kind)(t)
        qdes = act_matrix.T @ xdes[:nq]
        vdes = act_matrix.T @ xdes[nq:]
        
        q = act_matrix.T @ x[:nq]
        v = act_matrix.T @ x[nq:]
        
        e = qdes - q
        edot = vdes - v
        
        u_feedback = Kp @ e + Kd @ edot
        
        return u_openloop + u_feedback, e, edot

    def dynamics(t, x):
        q = x[:nq]
        v = x[nq:]
        
        if ts is None:
            tau = np.zeros((nv))
        else:
            u, _, _ = control(t, x)
            tau = act_matrix @ u
        
        a = pin.constraintDynamics(
            model, data_sim, q, v, tau, constraint_model, constraint_data, prox_settings
        )
        
        return np.concatenate([v, a])
    
    sol = solve_ivp(dynamics, 
                    [ts_sim[0], ts_sim[-1]],
                    x0, 
                    method=method,
                    t_eval=ts_sim)
    
    controls = []
    position_errors = []
    velocity_errors = []
    
    if ts is not None:
        for i in range(len(ts_sim)):
            t = ts_sim[i]
            x = sol.y[:,i]
            u, e, edot = control(t, x)
            
            controls.append(u)
            position_errors.append(e)
            velocity_errors.append(edot)
            
    position = sol.y[:nq,:].T
    velocity = sol.y[nq:,:].T
    
    return position, velocity, controls, position_errors, velocity_errors

if __name__ == "__main__":
    urdf_filename = "../../../Robots/digit-v3/digit-v3-armfixedspecific-floatingbase-springfixed.urdf"
    model = pin.buildModelFromUrdf(urdf_filename)
    data = model.createData()
    
    nq = model.nq
    nv = model.nv
    nu = 12

    act_matrix = np.zeros((nv, nu))
    act_matrix[[6,7,8,11,13,16, 21,22,23,26,28,31], :] = np.eye(nu)
                # left leg      right leg
    
    # Create dynamics
    constraint_models = []

    # left toe A closed loop
    pl1 = pin.SE3.Identity()
    pl1.translation = np.array([0.17 * 2, 0, 0])
    pl2 = pin.SE3.Identity()
    pl2.translation = np.array([0.0179, -0.009551, -0.054164])
    contact_model_ltA = pin.RigidConstraintModel(
        pin.ContactType.CONTACT_3D,   
        model,
        model.getJointId('left_A2'),
        pl1,
        model.getJointId('left_toe_roll'),
        pl2,
    )
    contact_model_ltA.corrector.Kp[:] = (100, 100, 100)
    contact_model_ltA.corrector.Kd[:] = (10, 10, 10)
    constraint_models.extend([contact_model_ltA])

    # left toe B closed loop
    pl1 = pin.SE3.Identity()
    pl1.translation = np.array([0.144 * 2, 0, 0])
    pl2 = pin.SE3.Identity()
    pl2.translation = np.array([-0.0181, -0.009551, -0.054164])
    contact_model_ltB = pin.RigidConstraintModel(
        pin.ContactType.CONTACT_3D,   
        model,
        model.getJointId('left_B2'),
        pl1,
        model.getJointId('left_toe_roll'),
        pl2,
    )
    contact_model_ltB.corrector.Kp[:] = (100, 100, 100)
    contact_model_ltB.corrector.Kd[:] = (10, 10, 10)
    constraint_models.extend([contact_model_ltB])

    # left knee-tarsus closed loop
    pl1 = pin.SE3.Identity()
    pl1.translation = np.array([0.25 * 2, 0, 0])

    # heel spring transformation
    pl2_1 = pin.SE3.Identity()
    pl2_1.rotation = pin.rpy.rpyToMatrix(np.array([math.radians(4.47), math.radians(0.32), math.radians(155.8)]))
    pl2_1.translation = np.array([-0.01766, -0.029456, 0.00104])

    pl2_2 = pin.SE3.Identity()
    pl2_2.translation = np.array([0.113789, -0.011056, 0])

    pl2 = pl2_1 * pl2_2
    contact_model_lkt = pin.RigidConstraintModel(
        pin.ContactType.CONTACT_3D,   
        model,
        model.getJointId('left_ach2'),
        pl1,
        model.getJointId('left_tarsus'),
        pl2,
    )
    contact_model_lkt.corrector.Kp[:] = (100, 100, 100)
    contact_model_lkt.corrector.Kd[:] = (10, 10, 10)
    constraint_models.extend([contact_model_lkt])

    # right toe A closed loop
    pl1 = pin.SE3.Identity()
    pl1.translation = np.array([0.17 * 2, 0, 0])
    pl2 = pin.SE3.Identity()
    pl2.translation = np.array([0.0179, 0.009551, -0.054164])
    contact_model_rtA = pin.RigidConstraintModel(
        pin.ContactType.CONTACT_3D,   
        model,
        model.getJointId('right_A2'),
        pl1,
        model.getJointId('right_toe_roll'),
        pl2,
    )
    contact_model_rtA.corrector.Kp[:] = (100, 100, 100)
    contact_model_rtA.corrector.Kd[:] = (10, 10, 10)
    constraint_models.extend([contact_model_rtA])

    # right toe B closed loop
    pl1 = pin.SE3.Identity()
    pl1.translation = np.array([0.144 * 2, 0, 0])
    pl2 = pin.SE3.Identity()
    pl2.translation = np.array([-0.0181, 0.009551, -0.054164])
    contact_model_rtB = pin.RigidConstraintModel(
        pin.ContactType.CONTACT_3D,   
        model,
        model.getJointId('right_B2'),
        pl1,
        model.getJointId('right_toe_roll'),
        pl2,
    )
    contact_model_rtB.corrector.Kp[:] = (100, 100, 100)
    contact_model_rtB.corrector.Kd[:] = (10, 10, 10)
    constraint_models.extend([contact_model_rtB])

    # right knee-tarsus closed loop
    pl1 = pin.SE3.Identity()
    pl1.translation = np.array([0.25 * 2, 0, 0])
    pl2_1 = pin.SE3.Identity()
    pl2_1.rotation = pin.rpy.rpyToMatrix(np.array([math.radians(-4.47), math.radians(0.32), math.radians(-155.8)]))
    pl2_1.translation = np.array([-0.01766, 0.029456, 0.00104])

    pl2_2 = pin.SE3.Identity()
    pl2_2.translation = np.array([0.113789, 0.011056, 0])

    pl2 = pl2_1 * pl2_2
    contact_model_rkt = pin.RigidConstraintModel(
        pin.ContactType.CONTACT_3D,   
        model,
        model.getJointId('right_ach2'),
        pl1,
        model.getJointId('right_tarsus'),
        pl2,
    )
    contact_model_rkt.corrector.Kp[:] = (100, 100, 100)
    contact_model_rkt.corrector.Kd[:] = (10, 10, 10)
    constraint_models.extend([contact_model_rkt])

    # left foot contact
    pl_leftfoot = pin.SE3.Identity()
    pl_leftfoot.rotation = np.array([[0, 1, 0], [-0.5, 0, np.sin(np.pi/3)], [np.sin(np.pi/3), 0, 0.5]])
    pl_leftfoot.translation = np.array([0, -0.05456, -0.0315])
    contact_model_lfc = pin.RigidConstraintModel(
        pin.ContactType.CONTACT_6D,   
        model,
        model.getJointId('left_toe_roll'),
        pl_leftfoot,
        0,
        data.oMi[model.getJointId('left_toe_roll')] * pl_leftfoot,
        pin.LOCAL_WORLD_ALIGNED
    )
    contact_model_lfc.corrector.Kp[:] = (0, 0, 100, 0, 0, 0)
    contact_model_lfc.corrector.Kd[:] = (50, 50, 50, 50, 50, 50)
    constraint_models.extend([contact_model_lfc])

    # right foot contact
    pl_rightfoot = pin.SE3.Identity()
    pl_rightfoot.rotation = np.array([[0, -1, 0], [0.5, 0, -np.sin(np.pi/3)], [np.sin(np.pi/3), 0, 0.5]])
    pl_rightfoot.translation = np.array([0, 0.05456, -0.0315])
    contact_model_rfc = pin.RigidConstraintModel(
        pin.ContactType.CONTACT_6D,   
        model,
        model.getJointId('right_toe_roll'),
        pl_rightfoot,
        0,
        data.oMi[model.getJointId('right_toe_roll')] * pl_rightfoot,
        pin.LOCAL_WORLD_ALIGNED
    )
    contact_model_rfc.corrector.Kp[:] = (0, 0, 100, 0, 0, 0)
    contact_model_rfc.corrector.Kd[:] = (50, 50, 50, 50, 50, 50)
    constraint_models.extend([contact_model_rfc])
    
    constraint_datas = [cm.createData() for cm in constraint_models]
        
    # forward simulation settings
    dt_sim = 5e-4
    T_ss = 0.35
        
    # load results from RAPTOR
    step_length = 0.8
    trajectories = np.loadtxt('../data/solution-digit-forward-' + str(step_length) + '.txt')

    ts_raptor = np.linspace(0, T_ss, len(trajectories))
    xs_raptor = np.zeros((len(ts_raptor), nq + nv)) 
    us_raptor = np.zeros((len(ts_raptor), nu))
    
    for i, states in enumerate(trajectories):
        q = states[:nv]
        v = states[nv:(2*nv)]
        u = states[(2*nv):]
        
        xs_raptor[i, :nq] = q
        xs_raptor[i, nq:] = v
        us_raptor[i, :] = u
        
    # initial conditions
    q0 = xs_raptor[0, :nq]
    v0 = xs_raptor[0, nq:]

    pin.forwardKinematics(model, data, q0)
    pin.updateFramePlacements(model, data)

    x0 = np.concatenate((q0, v0))
    u0 = np.zeros(nu)
        
    # note that we only focus on the first left support stage
    constraint_model_left = constraint_models[:-1]
    constraint_data_left = constraint_datas[:-1]
    
    ts_sim = np.arange(0, T_ss, dt_sim)
    
    Kp = np.diag(80.0 * np.ones(nu))
    # Kd = np.diag(0.05 * np.sqrt(80.0) * np.ones(nu))
    Kd = np.diag(5.0 * np.ones(nu))
    
    pos_sim, vel_sim, us_sim, e_sim, edot_sim = integrate(
        model, constraint_model_left, constraint_data_left,
        ts_sim, x0,
        ts_raptor, xs_raptor, us_raptor, act_matrix,
        Kp, Kd)
    
    # swing foot statistics
    model.addFrame(pin.Frame('right_foot', model.getJointId('right_toe_roll'), 0, pl_rightfoot, pin.FrameType.OP_FRAME))
    data = model.createData()
    
    RF_id = model.getFrameId("right_foot")
    pin.forwardKinematics(model, data, xs_raptor[-1][:nq])
    pin.updateFramePlacements(model, data)
    RF_placement = data.oMf[RF_id]
    step_length_opt = -RF_placement.translation[1]
    pin.forwardKinematics(model, data, pos_sim[-1])
    pin.updateFramePlacements(model, data)
    RF_placement = data.oMf[RF_id]
    step_length_sim = -RF_placement.translation[1]
    print(step_length_opt, step_length_sim)
    
    np.savetxt('../data/trajectory-digit-simulation.txt', pos_sim.T)
    
    scipy.io.savemat('../data/digit-simulation-' + str(step_length) + '.mat', 
                     {'ts_sim': ts_sim, 
                      'pos_sim': pos_sim, 
                      'vel_sim': vel_sim, 
                      'us_sim': us_sim, 
                      'e_sim': e_sim, 
                      'edot_sim': edot_sim,
                      'step_length_opt': step_length_opt,
                      'step_length_sim': step_length_sim})
    
    