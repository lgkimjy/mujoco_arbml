//
#include "ARBML.h"

CARBML::CARBML() :	id_body(0), id_body_parent(0), kinematic_chain(0), joint_type(0), joint_axis(0), joint_dir(0),
					mass_G(0.0), mass_lnk(0), rho_LCS(0), Rot_LCS2CoM(0), I_CoM_diag(0), pos0_offset(0), Rot0_offset(0),
					IsDynamicsDone(0)
{
	_del_T = 0;
	_g_const = 0;
	_actJnt_start_bodyID = -1;
	_FloatingBaseFlag = -1;

	e1.setZero();
	e2.setZero();
	e3.setZero();
	e1(0) = e2(1) = e3(2) = 1;

	p_B.setZero();
	quat_B.setZero();
	pdot_B.setZero();
	omega_B.setZero();
	varphi_B.setZero();

	quat_B(0) = 1.0;

	q.setZero();
	qdot.setZero();

	E3.setIdentity();
	R_B.setIdentity();
	T_B.setIdentity();
	I_actuator.setZero();
}



////////////////////////////////////////////////////////////////////////////////
//	Initialize Robot System Parameters
////////////////////////////////////////////////////////////////////////////////
void CARBML::initRobot(const mjModel* model)
{
	int error = 0;

	//////////	0. Check DoF & dimension of variables
	if (TOTAL_DOF == model->nv && NO_OF_BODY == model->nbody - 1) error = 0;


	/////	Set Flags
#ifdef _FLOATING_BASE
	if (ACTIVE_DOF == model->njnt - 1) error = 0;
	_FloatingBaseFlag = 1;		//	Floating-base body Flag
	_actJnt_start_bodyID = 1;	//	Active joint starting body ID
#else
	if (ACTIVE_DOF == model->njnt) error = 0;
	_FloatingBaseFlag = 0;		//	Fixed-base body Flag
	_actJnt_start_bodyID = 0;	//	Active joint starting body ID
#endif


	//	 Reload & reset capacity of vectors
	clearCapacity();

	if (error == 0) {
		//////////	1. Bring ALL data of robot from MuJoCo
		initRobotModel(model);
	}
	else {
		_ErrorMsg("Check DoF & Dimension of Global Var. in Config.h");
	}

	assignCapacity();

	cout << endl << "[InitRobotModel_MJCF] : Model Construction Complete !! " << endl;

#ifdef PRINT_MODEL_INFO
	outputSysInformation(model);
#endif
}


void CARBML::initRobotModel(const mjModel* mjmodel)
{
	int i, jointID;
	Eigen::Vector3d temp_vec;
	Eigen::Vector4d	temp_quat;


	_del_T = mjmodel->opt.timestep;			//	Sampling Time 
	_g_const = -mjmodel->opt.gravity[2];	//	Gravity

	//////////	Set const. body & joint parameters from XML files
	for (i = 0; i < NO_OF_BODY + 1; i++) {		//	= mjmodel->nbody
		jointID = mjmodel->body_jntadr[i];
		if (jointID > -1) {			//	-1 : No Joint !
			id_body.push_back(i - 1);
			id_body_parent.push_back(mjmodel->body_parentid[i] - 1);

			///	Joint type : free joint = 0, prismatic joint(slide) = 2, revolute joint(hinge) = 3 
			joint_type.push_back(mjmodel->jnt_type[jointID]);

			///	Const. local joint axis
			temp_vec = {(sysReal)mjmodel->jnt_axis[jointID * 3],
						(sysReal)mjmodel->jnt_axis[jointID * 3 + 1],
						(sysReal)mjmodel->jnt_axis[jointID * 3 + 2]};
			joint_axis.push_back(temp_vec);


			if (temp_vec == e1 || temp_vec == -e1) {
				joint_dir.push_back(X_AXIS);
			}
			else if (temp_vec == e2 || temp_vec == -e2) {
				joint_dir.push_back(Y_AXIS);
			}
			else if (temp_vec == e3 || temp_vec == -e3) {
				joint_dir.push_back(Z_AXIS);
			}


			///	Const. local position offset rel. to parent body
			temp_vec = {(sysReal)mjmodel->body_pos[i * 3],
						(sysReal)mjmodel->body_pos[i * 3 + 1],
						(sysReal)mjmodel->body_pos[i * 3 + 2]};
			pos0_offset.push_back(temp_vec);

			///	Const. orientation offset rel. to parent body
			temp_quat = {(sysReal)mjmodel->body_quat[i * 4],
						 (sysReal)mjmodel->body_quat[i * 4 + 1],
						 (sysReal)mjmodel->body_quat[i * 4 + 2],
						 (sysReal)mjmodel->body_quat[i * 4 + 3]};
			Rot0_offset.push_back(_Quat2Rot(temp_quat));

			///	Const. link/body mass
			mass_lnk.push_back((sysReal)mjmodel->body_mass[i]);

			///	Const. CoM position offset rel. to local link frame
			temp_vec = {(sysReal)mjmodel->body_ipos[i * 3],
						(sysReal)mjmodel->body_ipos[i * 3 + 1],
						(sysReal)mjmodel->body_ipos[i * 3 + 2]};
			rho_LCS.push_back(temp_vec);

			///	Const. rotation matrix btw link frame and principal axis of CoM
			temp_quat = {(sysReal)mjmodel->body_iquat[i * 4],
					 	 (sysReal)mjmodel->body_iquat[i * 4 + 1],
					     (sysReal)mjmodel->body_iquat[i * 4 + 2],
						 (sysReal)mjmodel->body_iquat[i * 4 + 3]};
			Rot_LCS2CoM.push_back(_Quat2Rot(temp_quat));

			///	Const. the principal moments of inertia @ CoM (diagonal !!)
			temp_vec = {(sysReal)mjmodel->body_inertia[i * 3],
						(sysReal)mjmodel->body_inertia[i * 3 + 1],
						(sysReal)mjmodel->body_inertia[i * 3 + 2]};
			I_CoM_diag.push_back(Eigen::DiagonalMatrix<double, 3>{temp_vec(0), temp_vec(1), temp_vec(2)});
		}
	}
	id_body.shrink_to_fit();
	id_body_parent.shrink_to_fit();
	joint_dir.shrink_to_fit();
	joint_axis.shrink_to_fit();
	joint_type.shrink_to_fit();
	pos0_offset.shrink_to_fit();
	Rot0_offset.shrink_to_fit();
	mass_lnk.shrink_to_fit();
	rho_LCS.shrink_to_fit();
	Rot_LCS2CoM.shrink_to_fit();
	I_CoM_diag.shrink_to_fit();


	//////////	Set body properties & initialize the elementary transformation matrix
	mass_G = 0;
	for (i = 0; i < id_body.size(); i++) {		//	id_body.size() = NO_OF_BODY
		body[i].setProperty(pos0_offset[i], Rot0_offset[i],							//	Const. local pos/ori. offset
							mass_lnk[i], rho_LCS[i], Rot_LCS2CoM[i], I_CoM_diag[i],	//	Dynamic parameters
							joint_type[i], joint_axis[i], joint_dir[i]);			//	Joint type & local axis
		mass_G += body[i].get_mass();
	}

	//////////	Set joint armature : Added inertia of actuator
	for (i = 0; i < ACTIVE_DOF; i++)	I_actuator(i) = mjmodel->dof_armature[DOF_BASEBODY + i];


	//////////	Set kinematic chain : kinematic_chain.size() = NO_OF_BODY
	int parent_ID = 0;
	vector<unsigned> temp_chain;

	for (i = 0; i < id_body.size(); i++) {	//	id_body.size() = NO_OF_BODY
		temp_chain.clear();
		temp_chain.insert(temp_chain.begin(), i);
		parent_ID = id_body_parent.at(i);

		while (parent_ID > _actJnt_start_bodyID - 1) {
			temp_chain.insert(temp_chain.begin(), parent_ID);
			parent_ID = id_body_parent[parent_ID];
		}
		temp_chain.shrink_to_fit();
		kinematic_chain.push_back(temp_chain);
	}
	kinematic_chain.shrink_to_fit();	//	kinematic_chain.size() = NO_OF_BODY


	//////////	Set Limit for Joint Range of Motion
	for (i = 0; i < mjmodel->njnt; i++) {
		if (mjmodel->jnt_limited[i] == 1) {
			id_limited_joint.push_back(i);
			q_min.push_back(mjmodel->jnt_range[i * 2]);
			q_max.push_back(mjmodel->jnt_range[i * 2 + 1]);
		}
	}
	q_max.shrink_to_fit();
	q_min.shrink_to_fit();
	id_limited_joint.shrink_to_fit();
}



void CARBML::outputSysInformation(const mjModel* mjmodel)
{
	int i;

	cout << endl << "================================================================================" << endl;

	cout.precision(3);
	cout << "Number of DoF (mjmodel->nv) = " << mjmodel->nv << endl;
	cout << "Number of Bodies (mjmodel->nbody) = " << mjmodel->nbody << endl << endl;

	cout << "##### CHECK 01 : NO_OF_BODY = mjmodel->njnt ?" << endl;
	cout << " + Number of Joints (mjmodel->njnt) = " << mjmodel->njnt << endl;
	cout << " + Number of Moving Bodies (NO_OF_BODY) = " << NO_OF_BODY << endl << endl;

	cout << "##### CHECK 02 : All data below = NO_OF_BODY ( " << NO_OF_BODY << " ) ?" << endl;
	cout << " + id_body.size() = " << id_body.size() << endl;
	cout << " + id_body_parent.size() = " << id_body_parent.size() << endl;
	cout << " + joint_type.size() = " << joint_type.size() << endl;
	cout << " + kinematic_chain.size() = " << kinematic_chain.size() << endl << endl;

	cout << "Number of Active Joints (ACTIVE_DOF) = " << ACTIVE_DOF << endl;
	cout << "Active Joint Starting Body ID = " << _actJnt_start_bodyID << endl << endl;

	cout << "Body ID : Parent Body ID : Joint Type : Joint Axis Dir. : " << "Local Joint Axis Vector" << endl;
	for (i = 0; i < id_body.size(); i++) {
		cout << "   " << id_body.at(i) << "\t\t" << id_body_parent.at(i) << "\t\t";
		cout << joint_type.at(i) << "\t      " << joint_dir.at(i) << "\t\t ";
		const Eigen::IOFormat fmt(2, Eigen::DontAlignCols, ", ", "", "", "");
		cout << "     [ " << joint_axis[i].transpose().format(fmt) << " ]"<< endl;
	}

	//cout << endl << "Zero position offset in Local Frame" << endl;
	//for (i = 0; i < pos0_offset.size(); i++) {
	//	cout << "[" << i << "] : ";
	//	pos0_offset[i].print(2);
	//}

	//cout << endl << "Zero rotation offset in Local Frame" << endl;
	//for (i = 0; i < Rot0_offset.size(); i++) {
	//	cout << "[" << i << "] : " << endl;
	//	Rot0_offset[i].print();
	//}

	//cout << endl << "Link Mass = " << mass_G << endl;
	//for (i = 0; i < mass_lnk.size(); i++)
	//	cout << "[" << i << "] : " << mass_lnk[i] << endl;

	//cout << endl << "Local Position of CoM" << endl;
	//for (i = 0; i < rho_LCS.size(); i++) {
	//	cout << "[" << i << "] : ";
	//	rho_LCS[i].print(2);
	//}

	//cout << endl << "Principal Link Inertia@CoM (Diagonal)" << endl;
	//for (i = 0; i < I_CoM_diag.size(); i++) {
	//	cout << "[" << i << "] : " << endl;
	//	I_CoM_diag[i].print();
	//}

	//cout << endl << "Rot_CoM_LCS" << endl;
	//for (i = 0; i < Rot_LCS2CoM.size(); i++) {
	//	cout << "[" << i << "] : " << endl;
	//	Rot_LCS2CoM[i].print();
	//}

	cout << endl << "================================================================================" << endl;

	cout << endl << "Kinematic Chain" << endl;
	for (i = 0; i < kinematic_chain.size(); i++) {
		cout << "Size of " << i << "-th Chain = " << kinematic_chain[i].size() << " : ";
		for (auto& it : kinematic_chain[i])
			cout << it << " ";
		cout << endl;
	}

	cout << endl << "================================================================================" << endl;

	cout.precision(5);
	cout << endl << "Number of Limited Joints = " << id_limited_joint.size() << endl;
	for (i = 0; i < id_limited_joint.size(); i++) {
		cout << "Limited Joint ID = " << id_limited_joint.at(i) << " : ";
		cout << "[" << fixed << q_min.at(i) << "\t" << q_max.at(i) << "]" << endl;
	}

	cout << endl << "================================================================================" << endl;
}

////////////////////	End of Initialization	////////////////////



////////////////////////////////////////////////////////////////////////////////
//	Compute Kinematic Core for other use !!
//		* Pose(position/rotation), Jacobian and its time derivative w.r.t {B}
//		* Must be CALLED BEFORE all others
////////////////////////////////////////////////////////////////////////////////
void CARBML::computeMotionCore()
{
	IsDynamicsDone = 0;

	Sk_varphi_B = Skew(varphi_B);

	/////	Compute	link position & rotation matrix, link CoM position w.r.t {B}
	getPose_BCS();

	/////	Compute linear, angular and link CoM Jacobian w.r.t {B}
	getJacob_BCS();

	/////	Compute time derivative of linear / angular / link CoM Jacobian w.r.t {B}
	getJacobDeriv_BCS();


	for (int i = 0; i < NO_OF_BODY; i++) {
		/////	1. For fixed-base body system : etadot_k = rvel_k, varphi_k = omega_b2lnk_BCS
		etadot_lnkCoM[i] = rvel_lnkCoM[i];
		varphi_lnk[i] = omega_b2lnk_BCS[i];

#ifdef _FLOATING_BASE
		/////	2. For floating-base body system :
		/////		* etadot_k = rvel_k + vsrphi_b x rpos_k, varphi_k = omega_b2lnk_BCS + varphi_b
		etadot_lnkCoM[i] += varphi_B.cross(rpos_lnkCoM[i]);
		varphi_lnk[i] += varphi_B;
#endif
	}
}



/////	Compute position & roation of all link frames w.r.t {B}
void CARBML::getPose_BCS()
{
#ifdef _FLOATING_BASE
	rpos_lnkCoM[0] = rho_LCS[0];
	Sk_rpos_lnkCoM[0] = Skew(rpos_lnkCoM[0]);
#endif

	/////	Forward Kinematics from base body(floating-base or fixed-base) frame
	for (int i = _actJnt_start_bodyID; i < NO_OF_BODY; i++) {
		body[i].ElemetaryTransformation(q(i - _actJnt_start_bodyID));

		if (id_body_parent[i] == _actJnt_start_bodyID - 1) {
			Rot_B2Lnk[i] = body[i].Rot_Link();
			rpos_lnk[i] = body[i].pos_Link();
			rpos_lnkCoM[i] = body[i].pos_CoM();
		}
		else {
			Rot_B2Lnk[i] = Rot_B2Lnk[id_body_parent[i]] * body[i].Rot_Link();
			rpos_lnk[i] = rpos_lnk[id_body_parent[i]] + Rot_B2Lnk[id_body_parent[i]] * body[i].pos_Link();
			rpos_lnkCoM[i] = rpos_lnk[id_body_parent[i]] + Rot_B2Lnk[id_body_parent[i]] * body[i].pos_CoM();
		}

		Sk_rpos_lnkCoM[i] = Skew(rpos_lnkCoM[i]);
	}
}



/////	Compute Jacobian(linear/angular) of all links w.r.t {B}
void CARBML::getJacob_BCS()
{
	int i, j, col;


	/////	Get joint motion axis w.r.t base body frame {B}
	for (i = 0; i < NO_OF_BODY; i++) {
		if ((body[i].JointAxisType() == RevoluteJoint) || (body[i].JointAxisType() == PrismaticJoint)) {
			//jntAxis_BCS[i] = Rot_B2Lnk[i].col(body[i].JointAxisDir());
			jntAxis_BCS[i] = Rot_B2Lnk[i] * joint_axis[i];
		}
		else {
			jntAxis_BCS[i].setZero();
		}
	}

	for (i = _actJnt_start_bodyID; i < NO_OF_BODY; i++) {
		/////	Compute Jacobian matrix (linear/angular & link CoM) w.r.t {B}
		for (unsigned& kk : kinematic_chain[i]) {
			col = kk - _actJnt_start_bodyID;
			if (body[kk].JointAxisType() == RevoluteJoint) {
				Jr_lnk_BCS[i].col(col) = jntAxis_BCS[kk];		//	Link rotational Jacobian

				tempVec3 = jntAxis_BCS[kk].cross(rpos_lnk[i] - rpos_lnk[kk]);
				Jp_lnk_BCS[i].col(col) = tempVec3;				//	Link linear Jacobian

				tempVec3 = jntAxis_BCS[kk].cross(rpos_lnkCoM[i] - rpos_lnk[kk]);
				J_lnkCoM_BCS[i].col(col) = tempVec3;				//	Link CoM Jacobian
			}
			else if (body[kk].JointAxisType() == PrismaticJoint) {
				Jp_lnk_BCS[i].col(col) = jntAxis_BCS[kk];		//	Link linear Jacobian (prismatic joint)
				J_lnkCoM_BCS[i].col(col) = jntAxis_BCS[kk];		//	Link CoM Jacobian (prismatic joint)
			}
		}

		/////	Compute linear/angular/CoM velocity of each link w.r.t {B}
		rvel_lnk[i].setZero();
		rvel_lnkCoM[i].setZero();
		omega_b2lnk_BCS[i].setZero();

		for (j = 0; j < DOF3; j++) {
			for (unsigned& kk : kinematic_chain[i]) {
				col = kk - _actJnt_start_bodyID;
				rvel_lnk[i](j) += Jp_lnk_BCS[i](j, col) * qdot(col);
				rvel_lnkCoM[i](j) += J_lnkCoM_BCS[i](j, col) * qdot(col);
				omega_b2lnk_BCS[i](j) += Jr_lnk_BCS[i](j, col) * qdot(col);
			}
		}
	}
}



/////	Compute time derivative of Jacobian(linear/angular) of all links w.r.t {B}
void CARBML::getJacobDeriv_BCS()
{
	int i, col;


	/////	Time derivative of motion axis w.r.t. {B} for Jdot !
	for (i = _actJnt_start_bodyID; i < NO_OF_BODY; i++) {
		jntAxisdot_BCS[i] = omega_b2lnk_BCS[i].cross(jntAxis_BCS[i]);
	}

	/////	Time derivative of link CoM Jacobian w.r.t. {B} expressed in {B}
	for (i = _actJnt_start_bodyID; i < NO_OF_BODY; i++) {
		/////	1. For fixed-base body system : Zdot_k = Jdot_BCS !
		for (unsigned& j : kinematic_chain[i]) {
			col = j - _actJnt_start_bodyID;
			if (body[j].JointAxisType() == RevoluteJoint) {
				Zdotr_lnk[i].col(col) = jntAxisdot_BCS[j];

				tempVec3 = jntAxisdot_BCS[j].cross(rpos_lnkCoM[i] - rpos_lnk[j])
							+ jntAxis_BCS[j].cross(rvel_lnkCoM[i] - rvel_lnk[j]);
				Zdot_lnkCoM[i].col(col) = tempVec3;
			}
			else if (body[j].JointAxisType() == PrismaticJoint) {
				Zdot_lnkCoM[i].col(col) = jntAxisdot_BCS[j];
			}
		}

#ifdef _FLOATING_BASE
		/////	2. For floating-base body system : Zdot_k = Jdot_BCS + varphi_B x J_BCS !
		for (int row = 0; row < DOF3; row++) {
			for (unsigned& kk : kinematic_chain[i]) {
				col = kk - _actJnt_start_bodyID;
				for (int l = 0; l < DOF3; l++) {
					Zdotr_lnk[i](row, col) += Sk_varphi_B(row, l) * Jr_lnk_BCS[i](l, col);
					Zdot_lnkCoM[i](row, col) += Sk_varphi_B(row, l) * J_lnkCoM_BCS[i](l, col);
				}
			}
		}
#endif
	}
}

////////////////////	End of Motion Core Functions	////////////////////



////////////////////////////////////////////////////////////////////////////////
//	Compute general end-effector position & rotation for A body w.r.t {I}
//		* pos_offset & Rot_offset are local position & rotation values !
//		* Core function -> link / link CoM / end-effector pose can be found !!
////////////////////////////////////////////////////////////////////////////////
void CARBML::getBodyPose(const int& ID_body, const Eigen::Vector3d& pos_offset, const Eigen::Matrix3d& Rot_offset,	\
						Eigen::Vector3d& pos_body, Eigen::Matrix3d& Rot_body)
{
	int i, j, row, col;

	//	R_A = R_b * R_bk * R_kA
	//	p_A = p_b + R_b * (r_k + R_bk * d_A)
	for (row = 0; row < DOF3; row++) {
		pos_body(row) = p_B(row);
		for (col = 0; col < DOF3; col++) {
			Rot_body(row, col) = 0;
			pos_body(row) += R_B(row, col) * rpos_lnk[ID_body](col);
			for (i = 0; i < DOF3; i++) {
				pos_body(row) += R_B(row, col) * Rot_B2Lnk[ID_body](col, i) * pos_offset(i);
				for (j = 0; j < DOF3; j++) {
					Rot_body(row, col) += R_B(row, i) * Rot_B2Lnk[ID_body](i, j) * Rot_offset(j, col);
				}
			}
		}
	}
}



/////	Compute position & rotation of A link frame(ID_body) : R_offset = E3, p_offset = 0 !
void CARBML::getLinkPose(const int& ID_body, Eigen::Vector3d& pos_lnk, Eigen::Matrix3d& Rot_lnk)
{
	getBodyPose(ID_body, Eigen::Vector3d{0, 0, 0}, E3, pos_lnk, Rot_lnk);
}


/////	Compute position & rotation of A link CoM(ID_body) : R_offset = R_lcs2CoM, p_offset = rho !
void CARBML::getLinkCoMPose(const int& ID_body, Eigen::Vector3d& pos_lnkCoM, Eigen::Matrix3d& Rot_lnkCoM)
{
	getBodyPose(ID_body, rho_LCS[ID_body], Rot_LCS2CoM[ID_body], pos_lnkCoM, Rot_lnkCoM);
}



////////////////////////////////////////////////////////////////////////////////
//	Compute Jacobian of A rigid body w.r.t {I}
//		* Jacobian matrix of A body w.r.t {I} !
//		* point_I is global position vector
//		* Core function -> link / link CoM / end-effector
////////////////////////////////////////////////////////////////////////////////

//////////	01. Linear Jacobian Only : J = 3 x N !!
void CARBML::getBodyJacob(const int& bodyID, const Eigen::Vector3d& point_I, Eigen::Matrix<double, DOF3, TOTAL_DOF>& Jp_body)
{
	int k, col, row;

	/////	A point on body w.r.t {B}
	rpos_A = Eigen::Transpose(R_B) * (point_I - p_B);

	if (_FloatingBaseFlag != 1 || bodyID != 0) {
		/////	Compute A body Jacobian w.r.t. {B}
		for (unsigned& j : kinematic_chain[bodyID]) {
			col = j - _actJnt_start_bodyID;
			if (body[j].JointAxisType() == RevoluteJoint) {		///	Revolute Joint
				tempVec3 = jntAxis_BCS[j].cross(rpos_A - rpos_lnk[j]);
				Jp_BCS.col(col) = tempVec3;			//	Linear Jacobian
			}
			else {												///	Prismatic Joint
				Jp_BCS.col(col) = jntAxis_BCS[j];
			}
		}

		/////	Compute a rigid body Jacobian (Active Joint Part) w.r.t {I}
		Jp_body.setZero();
		for (row = 0; row < DOF3; row++) {
			for (unsigned& kk : kinematic_chain[bodyID]) {
				col = kk - _actJnt_start_bodyID;
				for (k = 0; k < DOF3; k++) {
					Jp_body(row, DOF_BASEBODY + col) += R_B(row, k) * Jp_BCS(k, col);
				}
			}
		}
	}

#ifdef _FLOATING_BASE
	tempMat3x3 = -Skew(R_B * rpos_A);

	for (row = 0; row < DOF3; row++) {
		Jp_body(row, row) = 1;
		for (col = 0; col < DOF3; col++) {
			Jp_body(row, DOF3 + col) = tempMat3x3(row, col);
		}
	}
#endif
}


//////////	02. Linear & Angular Jacobian : J = 6 x N !!
void CARBML::getBodyJacob(const int& bodyID, const Eigen::Vector3d& point_I,						\
						Eigen::Matrix<double, DOF3, TOTAL_DOF>& Jp_body, Eigen::Matrix<double, DOF3, TOTAL_DOF>& Jr_body)
{
	int k, col, row;

	/////	A point on body w.r.t {B}
	rpos_A = Eigen::Transpose(R_B) * (point_I - p_B);

	if (_FloatingBaseFlag != 1 || bodyID != 0) {
		/////	Compute A body Jacobian w.r.t. {B}
		for (unsigned& j : kinematic_chain[bodyID]) {
			col = j - _actJnt_start_bodyID;
			if (body[j].JointAxisType() == RevoluteJoint) {		///	Revolute Joint
				Jr_BCS.col(col) = jntAxis_BCS[j];				//	Angular Jacobian

				tempVec3 = jntAxis_BCS[j].cross(rpos_A - rpos_lnk[j]);
				Jp_BCS.col(col) = tempVec3;						//	Linear Jacobian
			}
			else {												///	Prismatic Joint
				Jp_BCS.col(col) = jntAxis_BCS[j];
			}
		}

		/////	Compute a rigid body Jacobian (Active Joint Part) w.r.t {I}
		Jp_body.setZero();
		Jr_body.setZero();
		for (row = 0; row < DOF3; row++) {
			for (unsigned& kk : kinematic_chain[bodyID]) {
				col = kk - _actJnt_start_bodyID;
				for (k = 0; k < DOF3; k++) {
					Jp_body(row, DOF_BASEBODY + col) += R_B(row, k) * Jp_BCS(k, col);
					Jr_body(row, DOF_BASEBODY + col) += R_B(row, k) * Jr_BCS(k, col);
				}
			}
		}
	}

#ifdef _FLOATING_BASE
	tempMat3x3 = -Skew(R_B * rpos_A);

	for (row = 0; row < DOF3; row++) {
		Jp_body(row, row) = 1;
		Jr_body(row, row + DOF3) = 1;
		for (col = 0; col < DOF3; col++) {
			Jp_body(row, DOF3 + col) = tempMat3x3(row, col);
		}
	}
#endif
}



////////////////////////////////////////////////////////////////////////////////
//	Compute Jacobian time derivative of A rigid body w.r.t {I}
//		* Jacobian time derivative matrix of A body w.r.t {I} !
//		* point_I is global position vector
//	- getBodyJacob() is followed by getBodyJacobDeriv()
////////////////////////////////////////////////////////////////////////////////

//////////	01. Linear Part Only : Jdot = 3 x N !!
void CARBML::getBodyJacobDeriv(const int& bodyID, Eigen::Matrix<double, DOF3, TOTAL_DOF>& Jdotp_body)
{
	int j(0), k(0), l(0), col(0), row(0);

	rvel_A.setZero();

	if (_FloatingBaseFlag != 1 || bodyID != 0) {
		/////	Compute local body velocity w.r.t. {B} for Jdot_BCS
		for (j = 0; j < DOF3; j++) {
			for (unsigned& kk : kinematic_chain[bodyID]) {
				col = kk - _actJnt_start_bodyID;
				rvel_A(j) += Jp_BCS(j, col) * qdot(col);
			}
		}

		/////	Compute time derivative of local body Jacobian w.r.t {B} : Zdot_BCS
		for (unsigned& kk : kinematic_chain[bodyID]) {		//	For fixed-base : Zdot_BCS = Jdot_BCS !!
			col = kk - _actJnt_start_bodyID;
			if (body[kk].JointAxisType() == RevoluteJoint) {	//	Revolute Joint
				tempVec3 = jntAxisdot_BCS[kk].cross(rpos_A - rpos_lnk[kk])
						+ jntAxis_BCS[kk].cross(rvel_A - rvel_lnk[kk]);
				Zdotp_BCS.col(col) = tempVec3;
			}
			else {												//	Prismatic Joint
				Zdotp_BCS.col(col) = jntAxisdot_BCS[kk];
			}
		}

#ifdef _FLOATING_BASE		//	For floating-base : Zdot_BCS = varphi_B x J_BCS + Jdot_BCS !!
		for (row = 0; row < DOF3; row++) {
			for (unsigned& jj : kinematic_chain[bodyID]) {
				col = jj - _actJnt_start_bodyID;
				for (l = 0; l < DOF3; l++) {
					Zdotp_BCS(row, col) += Sk_varphi_B(row, l) * Jp_BCS(l, col);
				}
			}
		}
#endif

		/////	Compute time derivative of body Jacobian(Active Joint Part) w.r.t {I}
		Jdotp_body.setZero();
		for (row = 0; row < DOF3; row++) {
			for (unsigned& kk : kinematic_chain[bodyID]) {
				col = kk - _actJnt_start_bodyID;
				for (k = 0; k < DOF3; k++) {
					Jdotp_body(row, DOF_BASEBODY + col) += R_B(row, k) * Zdotp_BCS(k, col);
				}
			}
		}
	}

#ifdef _FLOATING_BASE
	tempVec3 = Sk_varphi_B * rpos_A + rvel_A;
	tempMat3x3 = -Skew(R_B * tempVec3);

	for (row = 0; row < DOF3; row++) {
		for (col = 0; col < DOF3; col++) {
			Jdotp_body(row, DOF3 + col) = tempMat3x3(row, col);
		}
	}
#endif
}


//////////	02. Linear & Angular Jacobian Derivative : Jdot = 6 x N !!
void CARBML::getBodyJacobDeriv(const int& bodyID, Eigen::Matrix<double, DOF3, TOTAL_DOF>& Jdotp_body, Eigen::Matrix<double, DOF3, TOTAL_DOF>& Jdotr_body)
{
	int k, col, row;

	rvel_A.setZero();

	if (_FloatingBaseFlag != 1 || bodyID != 0) {
		/////	Compute local body velocity w.r.t. {B} for Jdot_BCS
		for (int j = 0; j < DOF3; j++) {
			for (unsigned& kk : kinematic_chain[bodyID]) {
				col = kk - _actJnt_start_bodyID;
				rvel_A(j) += Jp_BCS(j, col) * qdot(col);
			}
		}

		/////	Compute time derivative of local body Jacobian w.r.t {B} : Zdot_BCS
		for (unsigned& j : kinematic_chain[bodyID]) {		//	For fixed-base : Zdot_BCS = Jdot_BCS !!
			col = j - _actJnt_start_bodyID;
			if (body[j].JointAxisType() == RevoluteJoint) {		//	Revolute Joint
				Zdotr_BCS.col(col) = jntAxisdot_BCS[j];

				tempVec3 = jntAxisdot_BCS[j].cross(rpos_A - rpos_lnk[j])
						+ jntAxis_BCS[j].cross(rvel_A - rvel_lnk[j]);
				Zdotp_BCS.col(col) = tempVec3;
			}
			else {												//	Prismatic Joint
				Zdotp_BCS.col(col) = jntAxisdot_BCS[j];
			}
		}

#ifdef _FLOATING_BASE		//	For floating-base : Zdot_BCS = varphi_B x J_BCS + Jdot_BCS !!
		for (row = 0; row < DOF3; row++) {
			for (unsigned& j : kinematic_chain[bodyID]) {
				col = j - _actJnt_start_bodyID;
				for (int l = 0; l < DOF3; l++) {
					Zdotp_BCS(row, col) += Sk_varphi_B(row, l) * Jp_BCS(l, col);
					Zdotr_BCS(row, col) += Sk_varphi_B(row, l) * Jr_BCS(l, col);
				}
			}
		}
#endif

		/////	Compute time derivative of body Jacobian(Active Joint Part) w.r.t {I}
		Jdotp_body.setZero();
		Jdotr_body.setZero();
		for (row = 0; row < DOF3; row++) {
			for (unsigned& kk : kinematic_chain[bodyID]) {
				col = kk - _actJnt_start_bodyID;
				for (k = 0; k < DOF3; k++) {
					Jdotp_body(row, DOF_BASEBODY + col) += R_B(row, k) * Zdotp_BCS(k, col);
					Jdotr_body(row, DOF_BASEBODY + col) += R_B(row, k) * Zdotr_BCS(k, col);
				}
			}
		}
	}

#ifdef _FLOATING_BASE
	tempVec3 = Sk_varphi_B * rpos_A + rvel_A;
	tempMat3x3 = -Skew(R_B * tempVec3);

	for (row = 0; row < DOF3; row++) {
		for (col = 0; col < DOF3; col++) {
			Jdotp_body(row, DOF3 + col) = tempMat3x3(row, col);
		}
	}
#endif
}



////////////////////////////////////////////////////////////////////////////////
//	Compute CoM kinematics of the system w.r.t {I}
//		* CoM position w.r.t {I}
//		* CoM Jacobian & its time derivative w.r.t {I} !
////////////////////////////////////////////////////////////////////////////////
void CARBML::computeCoMKinematics()
{
	getCoMPosition();		//	Compute p_B2CoM & p_CoM

	getCoMJacobian();		//	Compute CoM Jacobian matrix

	getCoMJacobianDeriv();

	pdot_CoM = J_CoM * xidot;
}


/////	Compute CoM position of the system
void CARBML::getCoMPosition()
{
	int i;

	/////	Compute CoM position of the system w.r.t {I}
	p_B2CoM = body[0].get_mass() * rpos_lnkCoM[0];
	for (i = 1; i < NO_OF_BODY; i++) {
		p_B2CoM += body[i].get_mass() * rpos_lnkCoM[i];
	}

	if (mass_G != 0)	p_B2CoM = p_B2CoM / mass_G;
	else				_ErrorMsg("Zero Mass !!!");

	p_B2CoM = R_B * p_B2CoM;
	p_CoM = p_B + p_B2CoM;
}


/////	Compute CoM Jacobian matrix of the system
void CARBML::getCoMJacobian()
{
	int i(0), j(0), k(0);
	Eigen::Matrix<double, DOF3, ACTIVE_DOF> J_CoM_BCS;

	J_CoM_BCS = body[0].get_mass() * J_lnkCoM_BCS[0];

	for (i = 1; i < NO_OF_BODY; i++) {
		for (j = 0; j < DOF3; j++) {
			for (unsigned& kk : kinematic_chain[i]) {
				k = kk - _actJnt_start_bodyID;
				J_CoM_BCS(j, k) += body[i].get_mass() * J_lnkCoM_BCS[i](j, k);
			}
		}
	}

	if (mass_G != 0)	J_CoM_BCS = J_CoM_BCS / mass_G;

	J_CoM.setZero();
	for (i = 0; i < DOF3; i++) {
		for (j = 0; j < ACTIVE_DOF; j++) {
			for (k = 0; k < DOF3; k++) {
				J_CoM(i, DOF_BASEBODY + j) += R_B(i, k) * J_CoM_BCS(k, j);
			}
		}
	}

#ifdef _FLOATING_BASE
	tempMat3x3 = -Skew(p_B2CoM);
	for (i = 0; i < DOF3; i++) {
		J_CoM(i, i) = 1;
		for (j = 0; j < DOF3; j++) {
			J_CoM(i, DOF3 + j) = tempMat3x3(i, j);
		}
	}
#endif
}


/////	Compute time derivative of CoM Jacobian matrix
void CARBML::getCoMJacobianDeriv()
{
	int i(0), j(0), k(0);
	Eigen::Vector3d etadot_G;
	Eigen::Matrix<double, DOF3, ACTIVE_DOF> Zdot_CoM_BCS;

	Zdot_CoM_BCS = body[0].get_mass() * Zdot_lnkCoM[0];
	etadot_G = body[0].get_mass() * etadot_lnkCoM[0];

	for (i = 1; i < NO_OF_BODY; i++) {
		for (j = 0; j < DOF3; j++) {
			etadot_G(j) += body[i].get_mass() * etadot_lnkCoM[i](j);
			for (unsigned& kk : kinematic_chain[i]) {
				k = kk - _actJnt_start_bodyID;
				Zdot_CoM_BCS(j, k) += body[i].get_mass() * Zdot_lnkCoM[i](j, k);
			}
		}
	}

	if (mass_G != 0) {
		pdot_B2CoM = R_B * etadot_G / mass_G;
		Zdot_CoM_BCS = Zdot_CoM_BCS / mass_G;
	}

	Jdot_CoM.setZero();
	for (i = 0; i < DOF3; i++) {
		for (j = 0; j < ACTIVE_DOF; j++) {
			for (k = 0; k < DOF3; k++) {
				Jdot_CoM(i, DOF_BASEBODY + j) += R_B(i, k) * Zdot_CoM_BCS(k, j);
			}
		}
	}

#ifdef _FLOATING_BASE
	tempMat3x3 = -Skew(pdot_B2CoM);
	for (i = 0; i < DOF3; i++) {
		for (j = 0; j < DOF3; j++) {
			Jdot_CoM(i, DOF3 + j) = tempMat3x3(i, j);
		}
	}
#endif
}



////////////////////////////////////////////////////////////////////////////////
//	Compute joint space dynamics
//		* xidot = ( pdot_b, omega_b, qdot )
////////////////////////////////////////////////////////////////////////////////
void CARBML::computeDynamics()
{
	IsDynamicsDone = 1;

	for (int i = 0; i < NO_OF_BODY; i++) {
		Sk_etadot_lnkCoM[i] = Skew(etadot_lnkCoM[i]);
		Sk_varphi_lnk[i] = Skew(varphi_lnk[i]);

		I_G_BCS[i] = Rot_B2Lnk[i] * body[i].Inertia_CoM() * Eigen::Transpose(Rot_B2Lnk[i]);
	}

#ifdef _FLOATING_BASE
	T_B.block(0, 0, R_B.rows(), R_B.cols()) = R_B;
	T_B.block(DOF3, DOF3, R_B.rows(), R_B.cols()) = R_B;
#endif

	getInertiaMatrix();
	getCoriolisCentrifugalMatrix();
	getGravityForce();

#ifdef _FLOATING_BASE
	computeCentroidalDynamics();
#endif

#if defined(INERTIADOT)
	getInertiaDot();
#endif
}


void CARBML::getInertiaMatrix()
{
	int i, j, k, l = 0, u = 0;
	int row = 0, col = 0;

	M_mat.setZero();

#ifdef _FLOATING_BASE	///	For floating-base body system : _actJnt_start_bodyID = 1
	for (i = 0; i < DOF3; i++) {
		M_mat(i, i) = mass_G;																	//	M11
	}

	/////	For all bodies !!
	for (i = 0; i < NO_OF_BODY; i++) {
		tempMat3x3 = body[i].get_mass() * Sk_rpos_lnkCoM[i] * Sk_rpos_lnkCoM[i];
		for (j = 0; j < DOF3; j++) {
			for (k = 0; k < DOF3; k++) {
				M_mat(j, DOF3 + k) -= body[i].get_mass() * Sk_rpos_lnkCoM[i](j, k);				//	M12
				M_mat(DOF3 + k, j) = M_mat(j, DOF3 + k);										//	M21 = M12^T
				M_mat(DOF3 + j, DOF3 + k) += (I_G_BCS[i](j, k) - tempMat3x3(j, k));				//	M22
			}
		}
	}

	/////	For active joints !!!
	for (i = _actJnt_start_bodyID; i < NO_OF_BODY; i++) {
		for (j = 0; j < DOF3; j++) {
			for (unsigned& kk : kinematic_chain[i]) {
				col = kk - _actJnt_start_bodyID;
				M_mat(j, DOF_BASEBODY + col) += body[i].get_mass() * J_lnkCoM_BCS[i](j, col);	//	M13
				M_mat(DOF_BASEBODY + col, j) = M_mat(j, DOF_BASEBODY + col);					//	M31 = M13^T
			}
		}

		tempMat3xAct_1.setZero();
		tempMat3xAct_2.setZero();
		for (j = 0; j < DOF3; j++) {
			for (unsigned& kk : kinematic_chain[i]) {
				col = kk - _actJnt_start_bodyID;
				for (l = 0; l < DOF3; l++) {
					tempMat3xAct_1(j, col) += I_G_BCS[i](j, l) * Jr_lnk_BCS[i](l, col);
					tempMat3xAct_2(j, col) += body[i].get_mass() * Sk_rpos_lnkCoM[i](j, l) * J_lnkCoM_BCS[i](l, col);
				}
				M_mat(DOF3 + j, DOF_BASEBODY + col) += (tempMat3xAct_1(j, col) + tempMat3xAct_2(j, col));	//	M_23
				M_mat(DOF_BASEBODY + col, DOF3 + j) = M_mat(DOF3 + j, DOF_BASEBODY + col);					//	M_32 = M23^T
			}
		}

		//	M_33
		for (unsigned& jj : kinematic_chain[i]) {
			row = jj - _actJnt_start_bodyID;
			for (unsigned& kk : kinematic_chain[i]) {
				col = kk - _actJnt_start_bodyID;
				for (l = 0; l < DOF3; l++) {
					for (u = 0; u < DOF3; u++) {
						M_mat(DOF_BASEBODY + row, DOF_BASEBODY + col) += Jr_lnk_BCS[i](l, row) * I_G_BCS[i](l, u) * Jr_lnk_BCS[i](u, col);
					}
					M_mat(DOF_BASEBODY + row, DOF_BASEBODY + col) += body[i].get_mass() * J_lnkCoM_BCS[i](l, row) * J_lnkCoM_BCS[i](l, col);
				}
			}
		}

		u = i - _actJnt_start_bodyID;
		M_mat(DOF_BASEBODY + u, DOF_BASEBODY + u) += I_actuator(u);
	}

	/////	Convert to the dynamics with absolute angular velocity of floating-base body
	M_mat = T_B * M_mat * Eigen::Transpose(T_B);

#else	///	For fixed-base body system : _actJnt_start_bodyID = 0
	for (i = 0; i < NO_OF_BODY; i++) {
		for (unsigned& row : kinematic_chain[i]) {
			for (unsigned& col : kinematic_chain[i]) {
				for (j = 0; j < DOF3; j++) {
					for (k = 0; k < DOF3; k++) {
						M_mat(row, col) += Jr_lnk_BCS[i](j, row) * I_G_BCS[i](j, k) * Jr_lnk_BCS[i](k, col);
					}
					M_mat(row, col) += body[i].get_mass() * J_lnkCoM_BCS[i](j, row) * J_lnkCoM_BCS[i](j, col);
				}
			}
		}

		/////	Add joint inertia of actuator : J_m * g_r * g_r
		M_mat(i, i) += I_actuator(i);
	}
#endif
}


void CARBML::getCoriolisCentrifugalMatrix()
{
	int i, j, k, l = 0, u = 0;
	int row = 0, col = 0;

	C_mat.setZero();

#ifdef _FLOATING_BASE	///	For floating-base body system : _actJnt_start_bodyID = 1
	/////	For floating-base body !
	Tmp_D = Sk_varphi_lnk[0] * I_G_BCS[0];
	tempMat3x3 = body[0].get_mass() * Sk_rpos_lnkCoM[0] * Sk_etadot_lnkCoM[0];

	C_mat.block(0, DOF3, DOF3, DOF3) = -body[0].get_mass() * Sk_etadot_lnkCoM[0];			//	C_12
	C_mat.block(DOF3, DOF3, DOF3, DOF3) = Tmp_D - tempMat3x3;								//	C_22

	/////	For active joint !
	for (i = _actJnt_start_bodyID; i < NO_OF_BODY; i++) {
		Tmp_D = Sk_varphi_lnk[i] * I_G_BCS[i];
		tempMat3x3 = body[i].get_mass() * Sk_rpos_lnkCoM[i] * Sk_etadot_lnkCoM[i];

		for (j = 0; j < DOF3; j++) {
			for (k = 0; k < DOF3; k++) {
				C_mat(j, DOF3 + k) -= body[i].get_mass() * Sk_etadot_lnkCoM[i](j, k);		//	C_12
				C_mat(DOF3 + j, DOF3 + k) += (Tmp_D(j, k) - tempMat3x3(j, k));				//	C_22
			}
		}

		for (j = 0; j < DOF3; j++) {
			for (unsigned& kk : kinematic_chain[i]) {
				col = kk - _actJnt_start_bodyID;
				C_mat(j, DOF_BASEBODY + col) += body[i].get_mass() * Zdot_lnkCoM[i](j, col);	//	C_13
			}
		}

		//	C_23
		for (j = 0; j < DOF3; j++) {
			for (unsigned& kk : kinematic_chain[i]) {
				col = kk - _actJnt_start_bodyID;
				for (l = 0; l < DOF3; l++) {
					C_mat(DOF3 + j, DOF_BASEBODY + col) += body[i].get_mass() * Sk_rpos_lnkCoM[i](j, l) * Zdot_lnkCoM[i](l, col)
														+ I_G_BCS[i](j, l) * Zdotr_lnk[i](l, col);
					C_mat(DOF3 + j, DOF_BASEBODY + col) += Tmp_D(j, l) * Jr_lnk_BCS[i](l, col);
				}
			}
		}

		//	C_32
		for (j = 0; j < DOF3; j++) {
			for (unsigned& kk : kinematic_chain[i]) {
				col = kk - _actJnt_start_bodyID;
				for (l = 0; l < DOF3; l++) {
					C_mat(DOF_BASEBODY + col, DOF3 + j) -= body[i].get_mass() * J_lnkCoM_BCS[i](l, col) * Sk_etadot_lnkCoM[i](l, j);
					C_mat(DOF_BASEBODY + col, DOF3 + j) += Jr_lnk_BCS[i](l, col) * Tmp_D(l, j);
				}
			}
		}

		//	C_33
		for (unsigned& jj : kinematic_chain[i]) {
			row = jj - _actJnt_start_bodyID;
			for (unsigned& kk : kinematic_chain[i]) {
				col = kk - _actJnt_start_bodyID;
				for (l = 0; l < DOF3; l++) {
					for (u = 0; u < DOF3; u++) {
						C_mat(DOF_BASEBODY + row, DOF_BASEBODY + col) += Jr_lnk_BCS[i](l, row) * (I_G_BCS[i](l, u) * Zdotr_lnk[i](u, col)
																	+ Tmp_D(l, u) * Jr_lnk_BCS[i](u, col));
					}
					C_mat(DOF_BASEBODY + row, DOF_BASEBODY + col) += body[i].get_mass() * J_lnkCoM_BCS[i](l, row) * Zdot_lnkCoM[i](l, col);
				}
			}
		}
	}

	/////	Convert to the dynamics with absolute angular velocity of floating-base body
	C_mat = T_B * C_mat * Eigen::Transpose(T_B);

#else	///	For fixed-base body system : _actJnt_start_bodyID = 0
	for (i = 0; i < NO_OF_BODY; i++) {
		Tmp_D = Sk_varphi_lnk[i] * I_G_BCS[i];
		for (unsigned& row : kinematic_chain[i]) {
			for (unsigned& col : kinematic_chain[i]) {
				for (j = 0; j < DOF3; j++) {
					for (k = 0; k < DOF3; k++) {
						C_mat(row, col) += Jr_lnk_BCS[i](j, row) * (I_G_BCS[i](j, k) * Zdotr_lnk[i](k, col)
										+ Tmp_D(j, k) * Jr_lnk_BCS[i](k, col));
					}
					C_mat(row, col) += body[i].get_mass() * J_lnkCoM_BCS[i](j, row) * Zdot_lnkCoM[i](j, col);
				}
			}
		}
	}
#endif
}


/////	Compute gravity force vector
void CARBML::getGravityForce()
{
	g_vec = mass_G * _g_const * J_CoM.row(Z_AXIS);
}



void CARBML::getInertiaDot()
{
	int i, l, u;
	int row = 0, col = 0;
	Eigen::Matrix3d Idot_G;	//	Time derivative of i-th link CoM Inertia matrix w.r.t. {B}
	Idot_G.setZero();		//	Time derivative of i-th link CoM Inertia matrix w.r.t. {B}

	Mdot_mat.setZero();

#ifdef _FLOATING_BASE	/////	Floating-base Body : _actJnt_start_bodyID = 1 !
	Idot_G = (Sk_varphi_lnk[0] * I_G_BCS[0] - I_G_BCS[0] * Sk_varphi_lnk[0]);
	tempMat3x3 = body[0].get_mass() * (Sk_etadot_lnkCoM[0] * Sk_rpos_lnkCoM[0] + Sk_rpos_lnkCoM[0] * Sk_etadot_lnkCoM[0]);

	for (int j = 0; j < DOF3; j++) {
		for (int k = 0; k < DOF3; k++) {
			Mdot_mat(j, DOF3 + k) -= body[0].get_mass() * Sk_etadot_lnkCoM[0](j, k);	//	Mdot12
			Mdot_mat(DOF3 + k, j) = Mdot_mat(j, DOF3 + k);								//	Mdot21 = Mdot12^T
			Mdot_mat(DOF3 + j, DOF3 + k) += (Idot_G(j, k) - tempMat3x3(j, k));			//	Mdot22
		}
	}

	/////	Active Joint !
	for (i = 1; i < NO_OF_BODY; i++) {
		Idot_G = (Sk_varphi_lnk[i] * I_G_BCS[i] - I_G_BCS[i] * Sk_varphi_lnk[i]);
		tempMat3x3 = body[i].get_mass() * (Sk_etadot_lnkCoM[i] * Sk_rpos_lnkCoM[i] + Sk_rpos_lnkCoM[i] * Sk_etadot_lnkCoM[i]);

		for (int j = 0; j < DOF3; j++) {
			for (int k = 0; k < DOF3; k++) {
				Mdot_mat(j, DOF3 + k) -= body[i].get_mass() * Sk_etadot_lnkCoM[i](j, k);	//	Mdot12
				Mdot_mat(DOF3 + k, j) = Mdot_mat(j, DOF3 + k);								//	Mdot21 = Mdot12^T
				Mdot_mat(DOF3 + j, DOF3 + k) += (Idot_G(j, k) - tempMat3x3(j, k));			//	Mdot22
			}
		}

		for (int j = 0; j < DOF3; j++) {
			for (unsigned& k : kinematic_chain[i]) {
				col = k - _actJnt_start_bodyID;
				Mdot_mat(j, DOF_BASEBODY + col) += body[i].get_mass() * Zdot_lnkCoM[i](j, col);	//	Mdot13
				Mdot_mat(DOF_BASEBODY + col, j) = Mdot_mat(j, DOF_BASEBODY + col);				//	Mdot31 = Mdot13^T
			}
		}

		tempMat3xAct_1.setZero();
		tempMat3xAct_2.setZero();
		for (int j = 0; j < DOF3; j++) {
			for (unsigned& k : kinematic_chain[i]) {
				col = k - _actJnt_start_bodyID;
				for (int l = 0; l < DOF3; l++) {
					tempMat3xAct_1(j, col) += Idot_G(j, l) * Jr_lnk_BCS[i](l, col) + I_G_BCS[i](j, l) * Zdotr_lnk[i](l, col);
					tempMat3xAct_2(j, col) += body[i].get_mass() * Sk_etadot_lnkCoM[i](j, l) * J_lnkCoM_BCS[i](l, col)
											+ body[i].get_mass() * Sk_rpos_lnkCoM[i](j, l) * Zdot_lnkCoM[i](l, col);
				}
				Mdot_mat(DOF3 + j, DOF_BASEBODY + col) += tempMat3xAct_1(j, col) + tempMat3xAct_2(j, col);	//	Mdot23
				Mdot_mat(DOF_BASEBODY + col, DOF3 + j) = Mdot_mat(DOF3 + j, DOF_BASEBODY + col);			//	Mdot32 = Mdot23^T
			}
		}

		//	Mdot33
		for (unsigned& j : kinematic_chain[i]) {
			row = j - _actJnt_start_bodyID;
			for (unsigned& k : kinematic_chain[i]) {
				col = k - _actJnt_start_bodyID;
				for (l = 0; l < DOF3; l++) {
					for (u = 0; u < DOF3; u++) {
						Mdot_mat(DOF_BASEBODY + row, DOF_BASEBODY + col) += Jr_lnk_BCS[i](l, row) * Idot_G(l, u) * Jr_lnk_BCS[i](u, col)
																		+ Zdotr_lnk[i](l, row) * I_G_BCS[i](l, u) * Jr_lnk_BCS[i](u, col)
																		+ Jr_lnk_BCS[i](l, row) * I_G_BCS[i](l, u) * Zdotr_lnk[i](u, col);
					}
					Mdot_mat(DOF_BASEBODY + row, DOF_BASEBODY + col) += body[i].get_mass() * Zdot_lnkCoM[i](l, row) * J_lnkCoM_BCS[i](l, col)
																	+ body[i].get_mass() * J_lnkCoM_BCS[i](l, row) * Zdot_lnkCoM[i](l, col);
				}
			}
		}
	}

	Mdot_mat = T_B * Mdot_mat * Eigen::Transpose(T_B);

#else	/////	Fixed-base Body System : _actJnt_start_bodyID = 0
	for (i = 0; i < NO_OF_BODY; i++) {
		Idot_G = Sk_varphi_lnk[i] * I_G_BCS[i] - I_G_BCS[i] * Sk_varphi_lnk[i];
		for (unsigned& j : kinematic_chain[i]) {
			row = j;
			for (unsigned& k : kinematic_chain[i]) {
				col = k;
				for (l = 0; l < DOF3; l++) {
					for (u = 0; u < DOF3; u++) {
						Mdot_mat(row, col) += Jr_lnk_BCS[i](l, row) * Idot_G(l, u) * Jr_lnk_BCS[i](u, col)
											+ Zdotr_lnk[i](l, row) * I_G_BCS[i](l, u) * Jr_lnk_BCS[i](u, col)
											+ Jr_lnk_BCS[i](l, row) * I_G_BCS[i](l, u) * Zdotr_lnk[i](u, col);
					}
					Mdot_mat(row, col) += body[i].get_mass() * Zdot_lnkCoM[i](l, row) * J_lnkCoM_BCS[i](l, col)
										+ body[i].get_mass() * J_lnkCoM_BCS[i](l, row) * Zdot_lnkCoM[i](l, col);
				}
			}
		}
	}
#endif
}


/////	Compute Centroidal Dynamics
void CARBML::computeCentroidalDynamics()
{
	Eigen::Matrix3d Sk_pos_BG;	//	Skew symmetric matrix form of pos_BG

	p_B2CoM(0) = M_mat(5, 1) / mass_G;
	p_B2CoM(1) = M_mat(3, 2) / mass_G;
	p_B2CoM(2) = M_mat(4, 0) / mass_G;

	Sk_pos_BG = Skew(p_B2CoM);

	M_p = M_mat.block(0, 0, DOF3, TOTAL_DOF);
	M_o = M_mat.block(DOF3, 0, DOF3, TOTAL_DOF);

	C_p = C_mat.block(0, 0, DOF3, TOTAL_DOF);
	C_o = C_mat.block(DOF3, 0, DOF3, TOTAL_DOF);

	l_b = M_p * xidot;
	k_b = M_o * xidot;
	h_b.segment(0, DOF3) = l_b;
	h_b.segment(DOF3, DOF3) = k_b;

	A_CoM.block(0, 0, DOF3, TOTAL_DOF) = M_p;
	A_CoM.block(DOF3, 0, DOF3, TOTAL_DOF) = M_o - Sk_pos_BG * M_p;

	h_CoM = A_CoM * xidot;
	l_CoM = h_CoM.segment(0, DOF3);
	k_CoM = h_CoM.segment(DOF3, DOF3);

	Adot_CoM.block(0, 0, DOF3, TOTAL_DOF) = C_p;
	Adot_CoM.block(DOF3, 0, DOF3, TOTAL_DOF) = C_o - Sk_pos_BG * C_p;
	
	hdot_CoM = A_CoM * xiddot + Adot_CoM * xidot;
	ldot_CoM = hdot_CoM.segment(0, DOF3);
	kdot_CoM = hdot_CoM.segment(DOF3, DOF3);

	pos_ZMP(X_AXIS) = p_CoM(X_AXIS) - ((p_CoM(Z_AXIS) * ldot_CoM(X_AXIS)) / (mass_G * _g_const + ldot_CoM(Z_AXIS)))
					- (kdot_CoM(Y_AXIS) / (mass_G * _g_const + ldot_CoM(Z_AXIS)));
	pos_ZMP(Y_AXIS) = p_CoM(Y_AXIS) - ((p_CoM(Z_AXIS) * ldot_CoM(Y_AXIS)) / (mass_G * _g_const + ldot_CoM(Z_AXIS)))
					+ (kdot_CoM(X_AXIS) / (mass_G * _g_const + ldot_CoM(Z_AXIS)));
	pos_ZMP(Z_AXIS) = 0.0;
}



////////////////////////////////////////////////////////////////////////////////
void CARBML::clearCapacity()
{
	/////	AssignCapacity()	/////
	id_body.clear();
	joint_axis.clear();
	joint_type.clear();
	joint_dir.clear();
	id_body_parent.clear();
	kinematic_chain.clear();

	rho_LCS.clear();
	mass_lnk.clear();
	Rot_LCS2CoM.clear();
	I_CoM_diag.clear();
	pos0_offset.clear();
	Rot0_offset.clear();

	q_max.clear();
	q_min.clear();
	id_limited_joint.clear();
	//////////////////////////////

	//J_lnkCoM.clear();
	//Jdot_lnkCoM.clear();
}


void CARBML::assignCapacity()
{
	int i;

	//J_lnkCoM.reserve(NO_OF_BODY);
	//Jdot_lnkCoM.reserve(NO_OF_BODY);

	for (i = 0; i < NO_OF_BODY; i++) {
		jntAxis_BCS[i].setZero();
		jntAxisdot_BCS[i].setZero();

		rpos_lnk[i].setZero();
		rpos_lnkCoM[i].setZero();
		rvel_lnk[i].setZero();
		rvel_lnkCoM[i].setZero();

		I_G_BCS[i].setZero();
		Rot_B2Lnk[i].setIdentity();

		varphi_lnk[i].setZero();
		etadot_lnkCoM[i].setZero();
		omega_b2lnk_BCS[i].setZero();

		Jp_lnk_BCS[i].setZero();
		Jr_lnk_BCS[i].setZero();
		J_lnkCoM_BCS[i].setZero();
		Zdotr_lnk[i].setZero();
		Zdot_lnkCoM[i].setZero();

		Sk_varphi_lnk[i].setZero();
		Sk_rpos_lnkCoM[i].setZero();
		Sk_etadot_lnkCoM[i].setZero();

		//J_lnkCoM[i].setZero();
		//Jdot_lnkCoM[i].setZero();
	}
}
