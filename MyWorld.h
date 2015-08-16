#ifndef _MYWORLD_
#define _MYWORLD_

#include <vector>
#include <Eigen/Dense>

class Particle;

class MyWorld {
 public:
    MyWorld(int _numParticles);

    virtual ~MyWorld();

    int getNumParticles() {
        return mParticles.size();
    }

    Particle* getParticle(int _index) {
        return mParticles[_index];
    }


    // TODO: your simulation code goes here
    void simulate(int mDisplayTimeout);
    Eigen::VectorXd getAccelerate();
    void setUpConstraints();
    void updateState();
    void midpoint(double h);
    void RK4(double h);

 protected:
    std::vector<Particle*> mParticles;
    int numOfConstraints;

    Eigen::MatrixXd J;
    Eigen::MatrixXd Jdot;
    Eigen::MatrixXd W;
    
    Eigen::VectorXd gravity;
    Eigen::VectorXd lambda;
    Eigen::VectorXd velocity;
    Eigen::VectorXd position;
    Eigen::VectorXd C;
    Eigen::VectorXd Cdot;
};

#endif
