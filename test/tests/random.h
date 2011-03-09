dmatrix m1(4, 4);
dvector v1(4);

mlcpp::UpdateRandomSeed();

std::cout << "mlcpp::RandVector<double>(5, 1) = " << mlcpp::RandVector<double>(5, 1) << std::endl;
std::cout << "mlcpp::RandVector<double>(5) = " << mlcpp::RandVector<double>(5) << std::endl;
std::cout << "mlcpp::RandMatrix<double>(5, 5, 1) = " << mlcpp::RandMatrix<double>(5, 5, 1) << std::endl;
std::cout << "mlcpp::RandMatrix<double>(5, 5) = " <<mlcpp::RandMatrix<double>(5, 5) << std::endl;

std::cout << "Randomize(m1, 10) = \n" << Randomize(m1, 10) << std::endl;
std::cout << "Randomize(m1) = \n" << Randomize(m1) << std::endl;
std::cout << "Randomize(v1, 10) = \n" << Randomize(v1, 10) << std::endl;
std::cout << "Randomize(v1) = \n" << Randomize(v1) << std::endl;

//std::cout << RandMatrix(5, 5, 1) << std::endl;
