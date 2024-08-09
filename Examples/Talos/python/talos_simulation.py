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
    urdf_filename = "../../../Robots/talos/talos_reduced_armfixed_floatingbase.urdf"
    model = pin.buildModelFromUrdf(urdf_filename)
    data = model.createData()
    
    nq = model.nq
    nv = model.nv
    nu = nv - 6
    
    FOOT_FRAME_IDS = {
        fname: model.getFrameId(fname) for fname in ["left_sole_link", "right_sole_link"]
    }
    FOOT_JOINT_IDS = {
        fname: model.frames[fid].parentJoint for fname, fid in FOOT_FRAME_IDS.items()
    }

    act_matrix = np.eye(nv, nu, -6)
    
    # Create dynamics
    prox_settings = pin.ProximalSettings(1e-9, 1e-10, 10)
    constraint_models = []
    constraint_datas = []
    for fname, fid in FOOT_FRAME_IDS.items():
        joint_id = FOOT_JOINT_IDS[fname]
        pl1 = model.frames[fid].placement
        pl2 = data.oMf[fid]
        cm = pin.RigidConstraintModel(
            pin.ContactType.CONTACT_6D,
            model,
            joint_id,
            pl1,
            0,
            pl2,
            pin.LOCAL_WORLD_ALIGNED,
        )
        cm.corrector.Kp[:] = (0, 0, 100, 0, 0, 0)
        cm.corrector.Kd[:] = (50, 50, 50, 50, 50, 50)
        constraint_models.append(cm)
        constraint_datas.append(cm.createData())
        
    # forward simulation settings
    dt_sim = 5e-4
    T_ss = 0.8
        
    # load results from RAPTOR
    step_length = 0.8
    trajectories = np.loadtxt('../data/solution-talos-forward-' + str(step_length) + '.txt')

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
    constraint_model_left = [constraint_models[0]]
    constraint_data_left = [constraint_datas[0]]
    
    ts_sim = np.arange(0, T_ss, dt_sim)
    
    Kp = np.diag(60.0 * np.ones(nu))
    # Kd = np.diag(0.05 * np.sqrt(60.0) * np.ones(nv))
    Kd = np.diag(5.0 * np.ones(nu))
    
    pos_sim, vel_sim, us_sim, e_sim, edot_sim = integrate(
        model, constraint_model_left, constraint_data_left,
        ts_sim, x0,
        ts_raptor, xs_raptor, us_raptor, act_matrix,
        Kp, Kd)
    
    RF_id = model.getFrameId("right_sole_link")
    pin.forwardKinematics(model, data, xs_raptor[-1][:nq])
    pin.updateFramePlacements(model, data)
    RF_placement = data.oMf[RF_id]
    step_length_opt = RF_placement.translation[0]
    pin.forwardKinematics(model, data, pos_sim[-1])
    pin.updateFramePlacements(model, data)
    RF_placement = data.oMf[RF_id]
    step_length_sim = RF_placement.translation[0]
    print(step_length_opt, step_length_sim)
    
    np.savetxt('../data/trajectory-talos-simulation.txt', pos_sim.T)
    
    scipy.io.savemat('../data/talos-simulation-' + str(step_length) + '.mat', 
                     {'ts_sim': ts_sim, 
                      'pos_sim': pos_sim, 
                      'vel_sim': vel_sim, 
                      'us_sim': us_sim, 
                      'e_sim': e_sim, 
                      'edot_sim': edot_sim,
                      'step_length_opt': step_length_opt,
                      'step_length_sim': step_length_sim})
    
    