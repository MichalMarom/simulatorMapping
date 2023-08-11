#include "WallDetection.h"

bool isNormallyDistributed(const std::vector<double>& data, double significance_level = 0.05) {
    if (data.empty()) 
    {
        // Empty data vector, cannot perform the test
        return false;
    }

    // Sort the data vector in ascending order
    std::vector<double> sorted_data = data;
    std::sort(sorted_data.begin(), sorted_data.end());

    // Calculate the mean and standard deviation of the data
    double mean = gsl_stats_mean(&sorted_data[0], 1, sorted_data.size());
    double stddev = gsl_stats_sd(&sorted_data[0], 1, sorted_data.size());

    // Calculate the test statistic (D) and p-value
    double D = 0.0;
    for (size_t i = 0; i < sorted_data.size(); ++i) 
    {
        double F_obs = gsl_cdf_ugaussian_P((sorted_data[i] - mean) / stddev);
        double F_exp = (i + 1.0) / sorted_data.size();
        double diff = std::abs(F_obs - F_exp);
        if (diff > D) 
        {
            D = diff;
        }
    }

    // Calculate the critical value for the given significance level and sample size
    double critical_value = gsl_cdf_ugaussian_Pinv(1.0 - significance_level / 2.0) / std::sqrt(sorted_data.size());

    return D <= critical_value;
}

/* Function that find the plane that minimizes the distance to a set of points
 input:
       vector<Eigen::Vector3d> points

 output:
        Eigen::Vector4d   plane equation (ax+by+cz+d=0)
*/
Eigen::Vector4d findMinimizingPlane(const vector<Eigen::Vector3d>& points) 
{
    Eigen::Matrix<double, Eigen::Dynamic, 3> A(points.size(), 3);

    // Fill the matrix A with points
    for (size_t i = 0; i < points.size(); ++i) 
    {
        A.row(i) = points[i].transpose();
    }

    // Compute the centroid of the points
    Eigen::Vector3d centroid = A.colwise().mean();

    // Subtract the centroid from each point to center them around the origin
    A.rowwise() -= centroid.transpose();

    // Compute the singular value decomposition (SVD) of A
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // The last column of V contains the normal vector of the plane
    Eigen::Vector3d normal = svd.matrixV().col(2);

    // The plane equation: ax + by + cz + d = 0
    // Compute the d value (distance from origin to the plane)
    double d = -normal.dot(centroid);

    // Return the plane as a 4D vector (a, b, c, d)
    return Eigen::Vector4d(normal.x(), normal.y(), normal.z(), d);
}

/* Function that calculate the angle between given plane to X - Z Plane
 input:
       Eigen::Vector3d& plane (coefficients a,b,c of the plane equation)

 output:
        double    degrees_angle (in degrees)   
*/
double angleBetweenPlanes(Eigen::Vector3d normalVector)
{
    // Normal vector of the X-Z plane
    Eigen::Vector3d xzPlaneNormal = Eigen::Vector3d(0, 1, 0); 

    double dot_prod = (normalVector[0] * xzPlaneNormal[0]) + (normalVector[1] * xzPlaneNormal[1]) + (normalVector[2] * xzPlaneNormal[2]);
    double mag1 = sqrt((normalVector[0] * normalVector[0]) + (normalVector[1] * normalVector[1]) + (normalVector[2] * normalVector[2]));
    double mag2 = sqrt((xzPlaneNormal[0] * xzPlaneNormal[0]) + (xzPlaneNormal[1] * xzPlaneNormal[1]) + (xzPlaneNormal[2] * xzPlaneNormal[2]));

    double cos_angle = dot_prod / (mag1 * mag2);
    double radian_angle = acos(cos_angle);

    // Convert the angle from radians to degrees
    double degrees_angle = radian_angle * (180.0 / M_PI);
    return degrees_angle;
}

/* Function that given points - decides whether they are a wall
 input:
       vector<Eigen::Vector3d> points

 output:
        bool  is_wall
*/
bool wallDetector(vector <Eigen::Vector3d>& points)
{
    bool is_wall = false;

    std::vector<double> z_cord;
    for (Eigen::Vector3d point : points) 
    {
        z_cord.push_back(point.z());
    }
    if (! isNormallyDistributed(z_cord))
    {
        is_wall = false;
        return is_wall;
    }

    // Find the plane that minimizes the distance to the points
    Eigen::Vector4d plane = findMinimizingPlane(points);
    Eigen::Vector3d plane_normal = Eigen::Vector3d(plane[0], plane[1], plane[2]);

    // Find the angle between the plane and XZ-plane
    double angle_between_plane_and_XZplane = angleBetweenPlanes(plane_normal);

    if (angle_between_plane_and_XZplane >= 88 && angle_between_plane_and_XZplane <= 92)
    {
        is_wall =  true;
    }
    else {
        is_wall = false;
    }
    return is_wall;
}
