#include "MyWorld.h"
#include "Particle.h"
#include <iostream>
using namespace Eigen;
using namespace std;

MyWorld::MyWorld(int _numParticles) {
  // Create particles
  for (int i = 0; i < _numParticles; i++) {
    Particle *p1 = new Particle();
    mParticles.push_back(p1);
   
  }
  cout<<"shuxin's tinkertoy"<<endl;
  // Init particle position
  mParticles[0]->mPosition[0] = 0.2;
  mParticles[1]->mPosition[0] = 0.2;
  mParticles[1]->mPosition[1] = -0.1;
  mParticles[2]->mPosition[0] = 0.2-0.05*sqrt(3);
  mParticles[2]->mPosition[1] = -0.05;

  numOfConstraints = 4;

  J.resize(numOfConstraints, mParticles.size()*3);
    
  Jdot.resize(numOfConstraints, mParticles.size()*3);
    
  W.resize(mParticles.size()*3, mParticles.size()*3);
  W.setIdentity();
  for(int i = 0; i < mParticles.size(); i++ )
		W.block(i*3, i*3, 3, 3) = Matrix3d::Identity() / mParticles[i]->mMass;

  
  gravity.resize(_numParticles*3, 1);
  for(int i = 0; i < _numParticles ; i++ )
     gravity.block<3,1>( i*3 , 0) << 0, -9.81, 0;
 
  velocity.resize(_numParticles*3, 1);
  for(int i = 0; i < _numParticles; i++ )
  	 velocity.block<3,1>( i*3, 0) << mParticles[i] -> mVelocity;

  position.resize(_numParticles*3, 1);
  for(int i = 0; i < _numParticles; i++ )
  	 position.block<3,1>( i*3, 0) << mParticles[i] -> mPosition;

  lambda.resize(_numParticles,1);
  C.resize(_numParticles,1);
  Cdot.resize(_numParticles, 1);
  setUpConstraints();
 } 

MyWorld::~MyWorld() {
  for (int i = 0; i < mParticles.size(); i++)
    delete mParticles[i];
  mParticles.clear();
}

void MyWorld::setUpConstraints(){
	C(0) = 0.5 * mParticles[0]->mPosition.dot(mParticles[0]->mPosition) - 0.5 * pow(0.2,2); 
  	C(1) = 0.5 * ( mParticles[0]->mPosition - mParticles[1]->mPosition).dot(mParticles[0]->mPosition - mParticles[1]->mPosition) - 0.5 * pow(0.1, 2);
  	C(2) = 0.5 * ( mParticles[1]->mPosition - mParticles[2]->mPosition).dot(mParticles[1]->mPosition - mParticles[2]->mPosition) - 0.5 * pow(0.1, 2);
  	C(3) = 0.5 * ( mParticles[0]->mPosition - mParticles[2]->mPosition).dot(mParticles[0]->mPosition - mParticles[2]->mPosition) - 0.5 * pow(0.1, 2);

  	Cdot(0) = mParticles[0] -> mPosition.dot (mParticles[0] -> mVelocity);
  	Cdot(1) = (mParticles[0] -> mPosition - mParticles[1] ->mPosition).dot(mParticles[0]->mVelocity - mParticles[1]->mVelocity);
  	Cdot(2) = (mParticles[1] -> mPosition - mParticles[2] ->mPosition).dot(mParticles[1]->mVelocity - mParticles[2]->mVelocity);
  	Cdot(3) = (mParticles[0] -> mPosition - mParticles[2] ->mPosition).dot(mParticles[0]->mVelocity - mParticles[2]->mVelocity);
}

VectorXd MyWorld::getAccelerate(){
	J.block(0,0,1,3) << position.block<3,1>(0,0).transpose();
	J.block(0,3,1,3) << 0 , 0, 0;
	J.block(0,6,1,3) << 0 , 0, 0;
	J.block(1,0,1,3) << (position.block<3,1>(0,0) - position.block<3,1>(3,0)).transpose();
	J.block(1,3,1,3) << (position.block<3,1>(3,0) - position.block<3,1>(0,0)).transpose();
	J.block(1,6,1,3) << 0 , 0, 0;
	J.block(2,0,1,3) << 0 , 0, 0;
	J.block(2,3,1,3) << (position.block<3,1>(3,0) - position.block<3,1>(6,0)).transpose();
	J.block(2,6,1,3) << (position.block<3,1>(6,0) - position.block<3,1>(3,0)).transpose();
	J.block(3,0,1,3) << (position.block<3,1>(0,0) - position.block<3,1>(6,0)).transpose();
	J.block(3,3,1,3) << 0 , 0, 0;
	J.block(3,6,1,3) << (position.block<3,1>(6,0) - position.block<3,1>(0,0)).transpose();
	

  

	Jdot.block(0,0,1,3) << velocity.block<3,1>(0,0).transpose();
	Jdot.block(0,3,1,3) << 0 , 0, 0;
	Jdot.block(0,6,1,3) << 0 , 0, 0;
	Jdot.block(1,0,1,3) << (velocity.block<3,1>(0,0) - velocity.block<3,1>(3,0)).transpose();
	Jdot.block(1,3,1,3) << (velocity.block<3,1>(3,0) - velocity.block<3,1>(0,0)).transpose();
	Jdot.block(1,6,1,3) << 0 , 0, 0;
	Jdot.block(2,0,1,3) << 0 , 0, 0;
	Jdot.block(2,3,1,3) << (velocity.block<3,1>(3,0) - velocity.block<3,1>(6,0)).transpose();
	Jdot.block(2,6,1,3) << (velocity.block<3,1>(6,0) - velocity.block<3,1>(3,0)).transpose();
	Jdot.block(3,0,1,3) << (velocity.block<3,1>(0,0) - velocity.block<3,1>(6,0)).transpose();
	Jdot.block(3,3,1,3) << 0 , 0, 0;
	Jdot.block(3,6,1,3) << (velocity.block<3,1>(6,0) - velocity.block<3,1>(0,0)).transpose();

 
	MatrixXd JWJT = J * W * J.transpose();
	MatrixXd B = -Jdot * velocity - J * W * gravity;
	lambda = JWJT.colPivHouseholderQr().solve(B);
	lambda += -4 * C - 4 * Cdot;
	VectorXd force = J.transpose() * lambda;
	
	return ( W * ( force + gravity) );//return local variable should return by copy

}
// integrators 
void MyWorld::RK4(double h){
	VectorXd initV = velocity;
	VectorXd initP = position;
	
	VectorXd A1 = getAccelerate();

	VectorXd P2 = initP + velocity * h * 0.5;
	VectorXd V2 = initV + A1 * h * 0.5;
	
	velocity = V2;
	position = P2;
	
	VectorXd A2 = getAccelerate();

	VectorXd P3 = initP + velocity * h * 0.5;
	VectorXd V3 = initV + A2 * h * 0.5;
	
	velocity = V3;
	position = P3;

	VectorXd A3 = getAccelerate();
	VectorXd P4 = initP+ velocity * h;
	VectorXd V4 = initV + A3 * h ;
	
	velocity = V4;
  	position = P4;

	VectorXd A4 = getAccelerate();

	VectorXd dxdt = 1.0/6.0 * initV + 1.0/3.0 * V2 + 1.0/3.0 * V3 + 1.0/6.0 * V4; 
	VectorXd dvdt = 1.0/6.0 * A1 + 1.0/3.0 * A2 + 1.0/3.0 * A3 + 1.0/6.0 * A4; 

	position = initP + dxdt * h;
	velocity = initV + dvdt * h;
}
void MyWorld::midpoint(double h){
	VectorXd initV = velocity;
	VectorXd initP = position;
	VectorXd A1 = getAccelerate();
	
	position = initP + velocity * h * 0.5;
	velocity = initV + A1 * h * 0.5;
	
	VectorXd A2 = getAccelerate();

	position = initP+ velocity * h ;
	velocity = initV + A2 * h ;
	
}
void MyWorld::updateState(){
	for(int i = 0; i < mParticles.size(); i++ )
  		mParticles[i] -> mVelocity << velocity.block<3,1>( i*3, 0);
  	
  	for(int i = 0; i < mParticles.size(); i++ )
  		mParticles[i] -> mPosition << position.block<3,1>( i*3, 0);
  	
}

void MyWorld::simulate(int mDisplayTimeout) {

	setUpConstraints();
	RK4(mDisplayTimeout / 1000.0);
	updateState();
	
}

