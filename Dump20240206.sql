-- MySQL dump 10.13  Distrib 8.0.32, for Linux (x86_64)
--
-- Host: localhost    Database: Quantum
-- ------------------------------------------------------
-- Server version	8.0.35-0ubuntu0.22.04.1

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!50503 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `cluster`
--

DROP TABLE IF EXISTS `cluster`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `cluster` (
  `idcluster` int NOT NULL,
  `grad` double DEFAULT NULL,
  `tugas` varchar(45) DEFAULT NULL,
  `ijk` varchar(45) DEFAULT NULL,
  `name` varchar(45) DEFAULT NULL,
  `grad_0` double DEFAULT NULL,
  PRIMARY KEY (`idcluster`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `cluster`
--

LOCK TABLES `cluster` WRITE;
/*!40000 ALTER TABLE `cluster` DISABLE KEYS */;
INSERT INTO `cluster` VALUES (0,-0.0006125955106580156,'1','0#0#0#0#','Na2,Na,Na}',-0.03751301711861288),(1,-0.0004946738082547299,'2','0#0#0#1#','Na2,Na,Na}',-0.03751301711861288),(2,-0.02547612641249778,'3','0#0#1#0#','Na2,Na,Na}',-0.03751301711861288),(3,-0.005149983280743973,'4','0#0#1#1#','Na2,Na,Na}',-0.03751301711861288),(4,-0.00021858365252280362,'5','0#0#2#0#','Na2,Na,Na}',-0.03751301711861288),(5,0.00016749871206959464,'1','0#0#2#1#','Na2,Na,Na}',-0.03751301711861288),(6,-0.018355453664149608,'2','0#1#0#0#','Na2,Na,Na}',-0.03751301711861288),(7,-0.021992782114925384,'3','0#1#0#1#','Na2,Na,Na}',-0.03751301711861288),(8,0.001923169755413674,'4','0#1#1#0#','Na2,Na,Na}',-0.03751301711861288),(9,-0.0064192945014119,'5','0#1#1#1#','Na2,Na,Na}',-0.03751301711861288),(10,-0.027229270817086336,'1','0#2#0#0#','Na2,Na,Na}',-0.03751301711861288),(11,-0.055706794032410516,'2','0#2#1#0#','Na2,Na,Na}',-0.03751301711861288),(12,-0.000012654882693509185,'3','1#0#0#0#','Na2,Na,Na}',-0.03751301711861288),(13,0.0000195346322625147,'4','1#0#0#1#','Na2,Na,Na}',-0.03751301711861288),(14,-0.05441963929354476,'5','1#0#1#0#','Na2,Na,Na}',-0.03751301711861288),(15,-0.05185970941789228,'1','1#0#1#1#','Na2,Na,Na}',-0.03751301711861288),(16,-0.000011640449327316303,'2','1#0#2#0#','Na2,Na,Na}',-0.03751301711861288),(17,0.00008317652139666631,'3','1#0#2#1#','Na2,Na,Na}',-0.03751301711861288),(18,-0.03767211385193325,'4','1#1#0#0#','Na2,Na,Na}',-0.03751301711861288),(19,-0.012737966887747606,'5','1#1#0#1#','Na2,Na,Na}',-0.03751301711861288),(20,NULL,'1','1#1#1#0#','Na2,Na,Na}',-0.03751301711861288),(21,NULL,'2','1#1#1#1#','Na2,Na,Na}',-0.03751301711861288),(22,NULL,'3','1#2#0#0#','Na2,Na,Na}',-0.03751301711861288),(23,NULL,'4','1#2#1#0#','Na2,Na,Na}',-0.03751301711861288),(24,NULL,'5','2#0#0#0#','Na2,Na,Na}',-0.03751301711861288),(25,NULL,'1','2#0#0#1#','Na2,Na,Na}',-0.03751301711861288),(26,NULL,'2','2#0#0#2#','Na2,Na,Na}',-0.03751301711861288),(27,NULL,'3','2#0#0#3#','Na2,Na,Na}',-0.03751301711861288);
/*!40000 ALTER TABLE `cluster` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `error`
--

DROP TABLE IF EXISTS `error`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `error` (
  `id` int NOT NULL AUTO_INCREMENT,
  `errorcol` double DEFAULT NULL,
  `from` varchar(45) DEFAULT NULL,
  `date` datetime DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=27 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `error`
--

LOCK TABLES `error` WRITE;
/*!40000 ALTER TABLE `error` DISABLE KEYS */;
INSERT INTO `error` VALUES (20,1.2227286922368277,NULL,'2024-02-06 14:48:28'),(21,1.2063665351924522,NULL,'2024-02-06 14:49:22'),(22,1.1330548712327269,NULL,'2024-02-06 14:50:28'),(23,1.8668832136700153,NULL,'2024-02-06 14:51:57'),(24,0.9017257725960899,NULL,'2024-02-06 14:53:07'),(25,1.0966800839688293,NULL,'2024-02-06 14:54:16'),(26,1.3554627733539066,NULL,'2024-02-06 14:55:23');
/*!40000 ALTER TABLE `error` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `mixer`
--

DROP TABLE IF EXISTS `mixer`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `mixer` (
  `Simbol` varchar(45) NOT NULL,
  `Energi_Ex` double DEFAULT NULL,
  `Energi_EJ` double DEFAULT NULL,
  `Energi_T` double DEFAULT NULL,
  `Energi_V` double DEFAULT NULL,
  `Energi_NRI` double DEFAULT NULL,
  PRIMARY KEY (`Simbol`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `mixer`
--

LOCK TABLES `mixer` WRITE;
/*!40000 ALTER TABLE `mixer` DISABLE KEYS */;
/*!40000 ALTER TABLE `mixer` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `mole`
--

DROP TABLE IF EXISTS `mole`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `mole` (
  `idmole` int NOT NULL AUTO_INCREMENT,
  `en` varchar(450) DEFAULT NULL,
  `molecol` varchar(450) DEFAULT NULL,
  PRIMARY KEY (`idmole`)
) ENGINE=InnoDB AUTO_INCREMENT=105 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `mole`
--

LOCK TABLES `mole` WRITE;
/*!40000 ALTER TABLE `mole` DISABLE KEYS */;
INSERT INTO `mole` VALUES (1,'0.0846165737051793','LiH (Lithium Hydride)'),(2,'0.0876736573705179','BeH (beryllium monohydride)'),(3,'0.131576573705179','CH (Methylidyne)'),(4,'0.488749992031873','CH3 (Methyl radical)'),(5,'0.660338007968127','CH4 (Methane)'),(6,'0.136278884462151','NH (Imidogen)'),(7,'0.290302422310757','NH2 (Amino radical)'),(8,'0.464789003984064','NH3 (Ammonia)'),(9,'0.167682709163347','OH (Hydroxyl radical)'),(10,'0.357569243027888','H2O (Water)'),(11,'0.215926533864542','HF (Hydrogen fluoride)'),(12,'0.23032812749004','SiH2 (silicon dihydride)'),(13,'0.349949800796813','SiH3 (Silyl radical)'),(14,'0.492521115537849','SiH4 (Silane)'),(15,'0.374370039840637','PH3 (Phosphine)'),(16,'0.283662151394422','H2S (Hydrogen sulfide)'),(17,'0.165988207171315','HCl (Hydrogen chloride)'),(18,'0.0328645418326693','Li2 (Lithium diatomic)'),(19,'0.212138326693227','LiF (lithium fluoride)'),(20,'0.646603075697211','C2H2 (Acetylene)'),(21,'0.895658358565737','C2H4 (Ethylene)'),(22,'1.12614565737052','C2H6 (Ethane)'),(23,'0.300099601593625','CN (Cyano radical)'),(24,'0.499361752988048','HCN (Hydrogen cyanide)'),(25,'0.408136892430279','CO (Carbon monoxide)'),(26,'0.452786135458167','HCO (Formyl radical)'),(27,'0.594453402390438','H2CO (Formaldehyde)'),(28,'0.805446215139442','CH3OH (Methyl alcohol)'),(29,'0.360158406374502','N2 (Nitrogen diatomic)'),(30,'0.685107888446215','N2H4 (Hydrazine)'),(31,'0.253950916334661','NO (Nitric oxide)'),(32,'0.216717131474104','O2 (Oxygen diatomic)'),(33,'0.417995219123506','H2O2 (Hydrogen peroxide)'),(34,'0.065123187250996','F2 (Fluorine diatomic)'),(35,'0.630339282868526','CO2 (Carbon dioxide)'),(36,'0.0294583266932271','Na2 (Disodium)'),(37,'0.126443187250996','Si2 (Silicon diatomic)'),(38,'0.183553147410359','P2 (Phosphorus diatomic)'),(39,'0.18049561752988','S2 (Sulfur diatomic)'),(40,'0.0958468525896414','Cl2 (Chlorine diatomic)'),(41,'0.149057370517928','NaCl (Sodium Chloride)'),(42,'0.292973067729084','SiO (Silicon monoxide)'),(43,'0.275004733067729','CS (carbon monosulfide)'),(44,'0.213678884462151','SO (Sulfur monoxide)'),(45,'0.119013545816733','ClO (Monochlorine monoxide)'),(46,'0.100738804780877','ClF (Chlorine monofluoride)'),(47,'0.816371314741036','Si2H6 (disilane)'),(48,'0.626029322709163','CH3Cl (Methyl chloride)'),(49,'0.260696812749004','HOCl (hypochlorous acid)'),(50,'0.408717928286853','SO2 (Sulfur dioxide)'),(51,')','BF3 (Borane trifluoro'),(52,')','BCl3 (Borane trichloro'),(53,'0.640957768924303','AlF3 (Aluminum trifluoride)'),(54,'0.476169402390438','AlCl3 (Aluminum trichloride)'),(55,'0.754655298804781','CF4 (Carbon tetrafluoride)'),(56,'0.514010358565737','CCl4 (Carbon tetrachloride)'),(57,'0.548621513944223','OCS (Carbonyl sulfide)'),(58,'0.464723505976096','CS2 (Carbon disulfide)'),(59,'0.67418809561753','CF2O (Carbonic difluoride)'),(60,'0.851760956175299','SiF4 (Silicon tetrafluoride)'),(61,')','SiCl4 (Silane tetrachloro'),(62,'0.455292430278884','N2O (Nitrous oxide)'),(63,'0.32717609561753','ClNO (Nitrosyl chloride)'),(64,'0.348519521912351','NF3 (Nitrogen trifluoride)'),(65,'0.256247011952191','O3 (Ozone)'),(66,'0.16446374501992','F2O (Difluorine monoxide)'),(67,'0.954019123505976','C2F4 (Tetrafluoroethylene)'),(68,')','CF3CN (Acetonitrile trifluoro'),(69,'1.12701721115538','CH3CCH (propyne)'),(70,'1.13262868525896','CH2CCH2 (allene)'),(71,'1.10061131474104','C3H4 (cyclopropene)'),(72,'1.36502039840637','C3H6 (Cyclopropane)'),(73,'1.59389752988048','C3H8 (Propane)'),(74,'1.59755697211155','C4H6 (Methylenecyclopropane)'),(75,'1.60722119521912','C4H6 (Cyclobutene)'),(76,'2.06172094023904','CH3CH(CH3)CH3 (Isobutane)'),(77,'2.20714103585657','C6H6 (Benzene)'),(78,')','CH2F2 (Methane difluoro'),(79,')','CHF3 (Methane trifluoro'),(80,'0.59290422310757','CH2Cl2 (Methylene chloride)'),(81,'0.55529561752988','CHCl3 (Chloroform)'),(82,'0.847816733067729','CH3CN (Acetonitrile)'),(83,'0.797303952191235','HCOOH (Formic acid)'),(84,'1.38274533864542','CH3CONH2 (Acetamide)'),(85,'1.15264637450199','C2H5N (Aziridine)'),(86,'0.818261354581673','C2N2 (Cyanogen)'),(87,'0.862510278884462','CH2CO (Ketene)'),(88,'1.04239378486056','C2H4O (Ethylene oxide)'),(89,'1.01647266932271','C2H2O2 (Ethanedial)'),(90,'1.27699505976096','CH3CH2OH (Ethanol)'),(91,'1.2594009561753','CH3OCH3 (Dimethyl ether)'),(92,'1.04239378486056','C2H4O (Ethylene oxide)'),(93,')','CH2CHF (Ethene fluoro'),(94,'1.0982035059761','CH3CH2Cl (Ethyl chloride)'),(95,'1.12709944223108','CH3COF (Acetyl fluoride)'),(96,'1.07315585657371','CH3COCl (Acetyl Chloride)'),(97,'1.60508207171315','C4H4O (Furan)'),(98,'1.73094183266932','C4H5N (Pyrrole)'),(99,'0.165396015936255','H2 (Hydrogen diatomic)'),(100,'0.137352828685259','HS (Mercapto radical)'),(101,'0.431237290836653','C2H (Ethynyl radical)'),(102,'0.641651314741036','CH3O (Methoxy radical)'),(103,'0.610272828685259','CH3S (thiomethoxy)'),(104,'0.396884462151394','NO2 (Nitrogen dioxide)');
/*!40000 ALTER TABLE `mole` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `pyscf`
--

DROP TABLE IF EXISTS `pyscf`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `pyscf` (
  `id` int NOT NULL AUTO_INCREMENT,
  `mole` varchar(45) DEFAULT NULL,
  `spin` varchar(45) DEFAULT NULL,
  `xyz` varchar(4500) DEFAULT NULL,
  `g2` varchar(45) DEFAULT NULL,
  `b3lyp` varchar(45) DEFAULT NULL,
  `blyp` varchar(45) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=63 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `pyscf`
--

LOCK TABLES `pyscf` WRITE;
/*!40000 ALTER TABLE `pyscf` DISABLE KEYS */;
INSERT INTO `pyscf` VALUES (2,'Al','1','13	0.0000	0.0000	0.0000','-241.930950','-242.365764','-242.350649'),(3,'Cl','1','17	0.0000	0.0000	0.0000','-459.676627','-460.131370','-460.112062'),(4,'H2','0','1	0.0000	0.0000	0.0000\n1	0.0000	0.0000	0.7414','-1.166358','-1.175482','-1.165255'),(5,'H','1','1	0.0000	0.0000	0.0000','-0.500000','-0.500273','-0.495446'),(6,'BeH','1','4	0.0000	0.0000	0.0000\n1	0.0000	0.0000	1.3426','-15.194853','-15.258221','-15.239766'),(7,'Be','0','4	0.0000	0.0000	0.0000','-14.622323','-14.668063','-14.655743'),(8,'LiH','0','3	0.0000	0.0000	0.0000\n1	0.0000	0.0000	1.5949','-8.022475','-8.080900','-8.065173'),(9,'Li','1','3	0.0000	0.0000	0.0000','-7.432217','-7.49090','-7.479999'),(10,'CH','1','6	0.0000	0.0000	0.0000\n1	0.0000	0.0000	1.1199','-38.412593','-38.469483','-38.451863'),(11,'C','2','6	0.0000	0.0000	0.0000','-37.78430','-37.843662','-37.829051'),(12,'CH3','1','6	0.0000	0.0000	0.0000\n1	1.0790	0.0000	0.0000\n1	-0.5395	-0.9344	0.0000\n1	-0.5395	0.9344	0.0000','-39.74508','-39.831174','-39.799139'),(13,'Cl2','0','17	0.0000	0.0000	0.0000\n17	0.0000	0.0000	1.9879','-919.442209','-920.317077','-920.288426'),(14,'H2O','0','8	0.0000	0.0000	0.1173\n1	0.0000	0.7572	-0.4692\n1	0.0000	-0.7572	-0.4692','-76.332049','-76.386117','-76.366676'),(15,'O','2','8	0.0000	0.0000	0.1173','-74.982030','-75.058330','-75.044562'),(16,'Li2','0','3	0.0000	0.0000	0.0000\n3	0.0000	0.0000	2.6730','-14.905760','-15.013884','-14.992012'),(17,'CH3OH','0','6	-0.0503	0.6685	0.0000\n8	-0.0503	-0.7585	0.0000\n1	-1.0807	1.0417	0.0000\n1	0.4650	1.0417	0.8924\n1	0.4650	1.0417	-0.8924\n1	0.8544	-1.0677	0.0000','-115.53489','-115.679502','-115.636418'),(18,'HCl','0','17	0.0000	0.0000	0.0000\n1	0.0000	0.0000	1.2746','-460.340177','-460.776281','-460.753551'),(19,'H2O2','0','8	0.0000	0.7375	-0.0528\n8	0.0000	-0.7375	-0.0528\n1	0.8190	0.8170	0.4220\n1	-0.8190	-0.8170	0.4220','-151.366513','-151.491088','-151.475124'),(20,'O2','2','8	0.0000	0.0000	0.0000\n8	0.0000	0.0000	1.2075','-150.148211','-150.270263','-150.272305'),(21,'CO','0','6	0.0000	0.0000	0.0000\n8	0.0000	0.0000	1.1282','-113.177498','-113.259179','-113.249096'),(22,'OH','1','8	0.0000	0.0000	0.0000\n1	0.0000	0.0000	0.9697','-75.64390','-75.708549','-75.693525'),(23,'CH4','0','6	0.0000	0.0000	0.0000\n1	0.6276	0.6276	0.6276\n1	0.6276	-0.6276	-0.6276\n1	-0.6276	0.6276	-0.6276\n1	-0.6276	-0.6276	0.6276','-40.410894','-40.510643','-40.472692'),(24,'NH','2','7	0.0000	0.0000	0.0000\n1	0.0000	0.0000	1.0362','-55.142178','-55.206969','-55.189577'),(25,'N','3','7	0.0000	0.0000	0.0000','-54.517960','-54.582875','-54.566371'),(26,'NH2','1','7	0.0000	0.0000	0.0000\n1	0.0000	0.8036	0.6347\n1	0.0000	-0.8036	0.6347','-55.78901','-55.852918','-55.831348'),(27,'NH3','0','7	0.0000	0.0000	0.0000\n1	0.0000	-0.9377	-0.3816\n1	0.8121	0.4689	-0.3816\n1	-0.8121	0.4689	-0.3816','-56.458638','-56.531887','-56.502486'),(28,'HF','0','9	0.0000	0.0000	0.0000\n1	0.0000	0.0000	0.9168','-100.350007','-100.403537','-100.388887'),(29,'F','1','9	0.0000	0.0000	0.0000','-99.632814','-99.713650','-99.700396'),(30,'LiF','0','9	0.0000	0.0000	0.0000\n3	0.0000	0.0000	1.5639','-107.284205','-107.404245','-107.38780'),(31,'C2H2','0','6	0.0000	0.0000	0.6013\n6	0.0000	0.0000	-0.6013\n1	0.0000	0.0000	1.6644\n1	0.0000	0.0000	-1.6644',NULL,'-77.31200','-77.279644'),(32,'C2H4','0','6	0.0000	0.0000	0.6695\n6	0.0000	0.0000	-0.6695\n1	0.0000	0.9289	1.2321\n1	0.0000	-0.9289	1.2321\n1	0.0000	0.9289	-1.2321\n1	0.0000	-0.9289	-1.2321','-78.415943','-78.572038','-78.523831'),(33,'C2H6','0','6	0.0000	0.0000	0.7680\n6	0.0000	0.0000	-0.7680\n1	-1.0192	0.0000	1.1573\n1	0.5096	0.8826	1.1573\n1	0.5096	-0.8826	1.1573\n1	1.0192	0.0000	-1.1573\n1	-0.5096	-0.8826	-1.1573\n1	-0.5096	0.8826	-1.1573','-79.631221','-79.812739','-79.74806'),(34,'CN','1','6	0.0000	0.0000	0.0000\n7	0.0000	0.0000	1.1718\n','-92.58275','-92.679753','-92.673648'),(35,'HCN','0','6	0.0000	0.0000	0\n7	0.0000	0.0000	1.064\n1	0.0000	0.0000	-1.156','-93.28489','-93.392489','-93.373276'),(36,'HCO','1','6	0.0000	0.0000	0.0000\n1	1.0800	0.0000	0.0000\n8	-0.5899	1.0427	0.0000','-113.698838','-113.807329','-113.792117'),(37,'H2CO','0','8	0.0000	0.0000	1.2050\n6	0.0000	0.0000	0.0000\n1	0.0000	0.9429	-0.5876\n1	0.0000	-0.9429	-0.5876\n','-114.338922','-114.461137','-114.437319'),(38,'N2','0','7	0.0000	0.0000	0.5488\n7	0.0000	0.0000	-0.5488',NULL,'-109.471040','-109.463192'),(39,'N2H4','0','7	0.0000	0.7230	-0.1123\n7	0.0000	-0.7230	-0.1123\n1	-0.4470	1.0031	0.7562\n1	0.4470	-1.0031	0.7562\n1	0.9663	1.0031	0.0301\n1	-0.9663	-1.0031	0.0301','-111.680423','-111.82001','-111.775885'),(40,'F2','0','9	0.0000	0.0000	0.0000\n9	0.0000	0.0000	1.4119\n','-199.323957','-199.477834','-199.476241'),(41,'CO2','0','6	0.0000	0.0000	0.0000\n8	0.0000	0.0000	1.1621\n8	0.0000	0.0000	-1.1621','-188.361321','-188.500268','-188.492599'),(42,'SiH3','1','14	0.0000	0.0000	0.0819\n1	0.0000	1.3928	-0.3820\n1	1.2062	-0.6964	-0.3820\n1	-1.2062	-0.6964	-0.3820\n','-290.773511','-291.206868','-291.173658'),(43,'Si','2','14	0.0000	0.0000	0.0819',NULL,'-289.368888','-289.352431'),(44,'SiH4','0','14	0.0000	0.0000	0.0000\n1	0.8544	0.8544	0.8544\n1	-0.8544	-0.8544	0.8544\n1	-0.8544	0.8544	-0.8544\n1	0.8544	-0.8544	-0.8544\n','-291.419051','-291.850760','-291.81127'),(45,'PH2','1','15	0.0000	0.0000	0.0000\n1	0.0000	1.0229	0.9964\n1	0.0000	-1.0229	0.9964\n','-342.049143','-342.477552','-342.451694'),(46,'P','3','15	0.0000	0.0000	0.0000','-340.818208','-341.255345','-341.236557'),(47,'PH3','0','15	0.0000	0.0000	0.0000\n1	0.0000	-1.1932	-0.7717\n1	1.0333	0.5966	-0.7717\n1	-1.0333	0.5966	-0.7717\n','-342.679023','-343.103832','-343.072607'),(48,'CH3SH','0','6	-0.8500	-0.0344	-0.2000\n16	0.9000	-0.5125	-0.1219\n1	1.4219	0.5781	0.4250\n1	-0.9406	0.8688	-0.8219\n1	-1.4219	-0.8688	-0.6469\n1	-1.2031	0.1656	0.8219\n','-438.148463','-438.660871','-438.608806'),(49,'S','2','16	0.9000	-0.5125	-0.1219','-397.65494','-398.100687','-398.08195'),(50,'H2S','0','16	0.0000	0.0000	0.1030\n1	0.0000	0.9616	-0.8239\n1	0.0000	-0.9616	-0.8239','-398.93071','-399.354995','-399.328984'),(51,'Na2','0','11	0.0000	0.0000	0.0000\n11	0.0000	0.0000	3.0789\n','-323.722985','-324.586697','-324.560023'),(52,'Na','1','11	0.0000	0.0000	0.0000','-161.846173','-162.279879','-162.266110'),(53,'Si2','2','14	0.0000	0.0000	0.0000\n14	0.0000	0.0000	2.2460\n',NULL,'-578.834257','-578.807115'),(54,'P2','0','15	0.0000	0.0000	0.0000\n15	0.0000	0.0000	1.8934\n','-681.819312','-682.642812','-682.623301'),(55,'S2','2','16	0.0000	0.0000	0.0000\n16	0.0000	0.0000	1.8892\n','-795.465114','-796.312703','-796.291691'),(56,'NaCl','0','11	0.0000	0.0000	0.0000\n17	0.0000	0.0000	2.3608\n','-621.680224','-622.550765','-622.516522'),(57,'CH3Cl','0','6	0.0000	0.0000	-1.1282\n17	0.0000	0.0000	0.6572\n1	0.0000	1.0357	-1.4679\n1	0.8969	-0.5179	-1.4679\n1	-0.8969	-0.5179	-1.4679\n','-499.553834','-500.081827','-500.033800'),(58,'FCl','0','9	0.0000	0.0000	0.0000\n17	0.0000	0.0000	1.6283','-559.406668','-559.912226','-559.895478'),(59,'SO2','0','16	0.0000	0.0000	0.0000\n8	0.0000	1.2371	0.7215\n8	0.0000	-1.2371	0.7215','-548.015738','-548.425599','-548.435036'),(60,'AlCl3','0','13	0.0000	0.0000	0.0000\n17	0.0000	2.0600	0.0000\n17	1.7840	-1.0300	0.0000\n17	-1.7840	-1.0300	0.0000','-1621.448747','-1623.143714','-1623.074280'),(61,'Al','1','13	0.0000	0.0000	0.0000','-241.930950','-242.365764','-242.350649'),(62,'SO','0','16	0.0000	0.0000	0.0000\n8	0.0000	0.0000	1.4811',NULL,NULL,NULL);
/*!40000 ALTER TABLE `pyscf` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `xyz_molecule`
--

DROP TABLE IF EXISTS `xyz_molecule`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `xyz_molecule` (
  `ID` int NOT NULL,
  `Nama` varchar(45) DEFAULT NULL,
  `Simbol` varchar(45) DEFAULT NULL,
  `Config` json DEFAULT NULL,
  `K_ref` double DEFAULT NULL,
  PRIMARY KEY (`ID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `xyz_molecule`
--

LOCK TABLES `xyz_molecule` WRITE;
/*!40000 ALTER TABLE `xyz_molecule` DISABLE KEYS */;
INSERT INTO `xyz_molecule` VALUES (1,'Beryllium monohydride','BeH','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 1.3426]}, \"atom\": [\"Be\", \"H\"], \"Spin_dn\": 2, \"Spin_up\": 3}',NULL),(2,'Lithium hydride','LiH','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 1.5949]}, \"atom\": [\"Li\", \"H\"], \"Spin_dn\": 2, \"Spin_up\": 2}',-8.022475),(3,'Hydrogen diatomic','H2','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 0.7414]}, \"atom\": [\"H\", \"H\"], \"Spin_dn\": 1, \"Spin_up\": 1}',-1.166358),(4,'Methylidyne','CH','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 1.1199]}, \"atom\": [\"C\", \"H\"], \"Spin_dn\": 3, \"Spin_up\": 4}',NULL),(5,'Methyl radical','CH3','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [1.079, 0, 0], \"3\": [-0.5395, -0.9344, 0], \"4\": [-0.5395, 0.9344, 0]}, \"atom\": [\"C\", \"H\", \"H\", \"H\"], \"Spin_dn\": 4, \"Spin_up\": 5}',-39.745085),(6,'Lithium','Li','{\"xyz\": {\"1\": [0, 0, 0]}, \"atom\": [\"Li\"], \"Spin_dn\": 2, \"Spin_up\": 1}',-7.432217),(7,'Beryllium','Be','{\"xyz\": {\"1\": [0, 0, 0]}, \"atom\": [\"Be\"], \"Spin_dn\": 2, \"Spin_up\": 2}',-14.622323),(8,'Chlorine diatomic','Cl2','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 1.9879]}, \"atom\": [\"Cl\", \"Cl\"], \"Spin_dn\": 17, \"Spin_up\": 17}',-919.442209),(9,'Hydrogen','H','{\"xyz\": {\"1\": [0, 0, 0]}, \"atom\": [\"H\"], \"Spin_dn\": 1, \"Spin_up\": 0}',-0.5),(10,'Helium','He','{\"xyz\": {\"1\": [0, 0, 0]}, \"atom\": [\"He\"], \"Spin_dn\": 1, \"Spin_up\": 1}',-2.900262),(11,'Carbon','C','{\"xyz\": {\"1\": [0, 0, 0]}, \"atom\": [\"C\"], \"Spin_dn\": 2, \"Spin_up\": 4}',NULL),(12,'Water','H2O','{\"xyz\": {\"1\": [0, 0, 0.1173], \"2\": [0, 0.7572, -0.4692], \"3\": [0, -0.7572, -0.4692]}, \"atom\": [\"O\", \"H\", \"H\"], \"Spin_dn\": 5, \"Spin_up\": 5}',NULL),(13,'Oksigen','O','{\"xyz\": {\"1\": [0, 0, 0]}, \"atom\": [\"O\"], \"Spin_dn\": 3, \"Spin_up\": 5}',NULL),(16,'Clorida','Cl','{\"xyz\": {\"1\": [0, 0, 0]}, \"atom\": [\"Cl\"], \"Spin_dn\": 9, \"Spin_up\": 8}',NULL),(18,'Litium','Li2','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 2.673]}, \"atom\": [\"Li\", \"Li\"], \"Spin_dn\": 3, \"Spin_up\": 3}',NULL),(19,'Methyl alcohol','CH3OH','{\"xyz\": {\"1\": [-0.0503, 0.6685, 0], \"2\": [-0.0503, -0.7585, 0], \"3\": [-1.0807, 1.0417, 0], \"4\": [0.465, 1.0417, 0.8924], \"5\": [0.465, 1.0417, -0.8924], \"6\": [0.8544, -1.0677, 0]}, \"atom\": [\"C\", \"O\", \"H\", \"H\", \"H\", \"H\"], \"Spin_dn\": 9, \"Spin_up\": 9}',NULL),(20,'Hydrogen chloride','HCl','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 1.2746]}, \"atom\": [\"H\", \"Cl\"], \"Spin_dn\": 9, \"Spin_up\": 9}',NULL),(22,'H2O2','H2O2','{\"xyz\": {\"1\": [0, 0.7375, -0.0528], \"2\": [0, -0.7375, -0.0528], \"3\": [0.819, 0.817, 0.422], \"4\": [-0.819, -0.817, 0.422]}, \"atom\": [\"O\", \"O\", \"H\", \"H\"], \"Spin_dn\": 9, \"Spin_up\": 9}',NULL),(23,'O2','O2','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 1.2075]}, \"atom\": [\"O\", \"O\"], \"Spin_dn\": 7, \"Spin_up\": 9}',NULL),(24,'CO','CO','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 1.1282]}, \"atom\": [\"C\", \"O\"], \"Spin_dn\": 7, \"Spin_up\": 7}',NULL),(25,'OH','OH','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 0.9697]}, \"atom\": [\"O\", \"H\"], \"Spin_dn\": 4, \"Spin_up\": 5}',NULL),(26,'CH4','CH4','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0.6276, 0.6276, 0.6276], \"3\": [0.6276, -0.6276, -0.6276], \"4\": [-0.6276, 0.6276, -0.6276], \"5\": [-0.6276, -0.6276, 0.6276]}, \"atom\": [\"C\", \"H\", \"H\", \"H\", \"H\"], \"Spin_dn\": 5, \"Spin_up\": 5}',NULL),(27,'NH','NH','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 1.0362]}, \"atom\": [\"N\", \"H\"], \"Spin_dn\": 3, \"Spin_up\": 5}',NULL),(28,'N','N','{\"xyz\": {\"1\": [0, 0, 0]}, \"atom\": [\"N\"], \"Spin_dn\": 2, \"Spin_up\": 5}',NULL),(30,'NH2','NH2','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0.8036, 0.6347], \"3\": [0, -0.8036, 0.6347]}, \"atom\": [\"N\", \"H\", \"H\"], \"Spin_dn\": 5, \"Spin_up\": 4}',NULL),(31,'NH3','NH3','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, -0.9377, -0.3816], \"3\": [0.8121, 0.4689, -0.3816], \"4\": [-0.8121, 0.4689, -0.3816]}, \"atom\": [\"N\", \"H\", \"H\", \"H\"], \"Spin_dn\": 5, \"Spin_up\": 5}',NULL),(32,'HF','HF','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 0.9168]}, \"atom\": [\"F\", \"H\"], \"Spin_dn\": 5, \"Spin_up\": 5}',NULL),(33,'F','F','{\"xyz\": {\"1\": [0, 0, 0]}, \"atom\": [\"F\"], \"Spin_dn\": 4, \"Spin_up\": 5}',NULL),(34,'LiF','LiF','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 1.5639]}, \"atom\": [\"F\", \"Li\"], \"Spin_dn\": 6, \"Spin_up\": 6}',NULL),(35,'C2H2','C2H2','{\"xyz\": {\"1\": [0, 0, 0.6013], \"2\": [0, 0, -0.6013], \"3\": [0, 0, 1.6644], \"4\": [0, 0, -1.6644]}, \"atom\": [\"C\", \"C\", \"H\", \"H\"], \"Spin_dn\": 7, \"Spin_up\": 7}',NULL),(36,'C2H4','C2H4','{\"xyz\": {\"1\": [0, 0, 0.6695], \"2\": [0, 0, -0.6695], \"3\": [0, 0.9289, 1.2321], \"4\": [0, -0.9289, 1.2321], \"5\": [0, 0.9289, -1.2321], \"6\": [0, -0.9289, -1.2321]}, \"atom\": [\"C\", \"C\", \"H\", \"H\", \"H\", \"H\"], \"Spin_dn\": 8, \"Spin_up\": 8}',NULL),(37,'C2H6','C2H6','{\"xyz\": {\"1\": [0, 0, 0.768], \"2\": [0, 0, -0.768], \"3\": [-1.0192, 0, 1.1573], \"4\": [0.5096, 0.8826, 1.1573], \"5\": [0.5096, -0.8826, 1.1573], \"6\": [1.0192, 0, -1.1573], \"7\": [-0.5096, -0.8826, -1.1573], \"8\": [-0.5096, 0.8826, -1.1573]}, \"atom\": [\"C\", \"C\", \"H\", \"H\", \"H\", \"H\", \"H\", \"H\"], \"Spin_dn\": 9, \"Spin_up\": 9}',NULL),(38,'CN','CN','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 1.1718]}, \"atom\": [\"C\", \"N\"], \"Spin_dn\": 6, \"Spin_up\": 7}',NULL),(39,'HCN','HCN','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 1.064], \"3\": [0, 0, -1.156]}, \"atom\": [\"C\", \"H\", \"N\"], \"Spin_dn\": 7, \"Spin_up\": 7}',NULL),(40,'HCO','HCO','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [1.08, 0, 0], \"3\": [-0.5899, 1.0427, 0]}, \"atom\": [\"C\", \"H\", \"O\"], \"Spin_dn\": 7, \"Spin_up\": 8}',NULL),(41,'H2CO','H2CO','{\"xyz\": {\"1\": [0, 0, 1.205], \"2\": [0, 0, 0], \"3\": [0, 0.9429, -0.5876], \"4\": [0, -0.9429, -0.587]}, \"atom\": [\"O\", \"C\", \"H\", \"H\"], \"Spin_dn\": 8, \"Spin_up\": 8}',NULL),(42,'N2','N2','{\"xyz\": {\"1\": [0, 0, -0.5488], \"2\": [0, 0, 0.5488]}, \"atom\": [\"N\", \"N\"], \"Spin_dn\": 7, \"Spin_up\": 7}',NULL),(43,'N2H4','N2H4','{\"xyz\": {\"1\": [0, 0.723, -0.1123], \"2\": [0, -0.723, -0.1123], \"3\": [-0.447, 1.0031, 0.7562], \"4\": [0.447, -1.0031, 0.7562], \"5\": [0.9663, 1.0031, 0.0301], \"6\": [-0.9663, -1.0031, 0.0301]}, \"atom\": [\"N\", \"N\", \"H\", \"H\", \"H\", \"H\"], \"Spin_dn\": 9, \"Spin_up\": 9}',NULL),(44,'F2','F2','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 1.4119]}, \"atom\": [\"F\", \"F\"], \"Spin_dn\": 9, \"Spin_up\": 9}',NULL),(45,'CO2','CO2','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 1.1621], \"3\": [0, 0, -1.1621]}, \"atom\": [\"C\", \"O\", \"O\"], \"Spin_dn\": 11, \"Spin_up\": 11}',NULL),(46,'SiH3','SiH3','{\"xyz\": {\"1\": [0, 0, 0.0819], \"2\": [0, 1.3928, -0.382], \"3\": [1.2062, -0.6964, -0.382], \"4\": [-1.2062, -0.6964, -0.382]}, \"atom\": [\"Si\", \"H\", \"H\", \"H\"], \"Spin_dn\": 8, \"Spin_up\": 9}',NULL),(47,'Si','Si','{\"xyz\": {\"1\": [0, 0, 0]}, \"atom\": [\"Si\"], \"Spin_dn\": 6, \"Spin_up\": 8}',NULL),(48,'SiH4','SiH4','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0.8544, 0.8544, 0.8544], \"3\": [-0.8544, -0.8544, 0.8544], \"4\": [-0.8544, 0.8544, -0.8544], \"5\": [0.8544, -0.8544, -0.8544]}, \"atom\": [\"Si\", \"H\", \"H\", \"H\", \"H\"], \"Spin_dn\": 9, \"Spin_up\": 9}',NULL),(49,'PH2','PH2','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 1.0229, 0.9964], \"3\": [0, -1.0229, 0.9964]}, \"atom\": [\"P\", \"H\", \"H\"], \"Spin_dn\": 8, \"Spin_up\": 9}',NULL),(50,'P','P','{\"xyz\": {\"1\": [0, 0, 0]}, \"atom\": [\"P\"], \"Spin_dn\": 6, \"Spin_up\": 9}',NULL),(51,'PH3','PH3','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, -1.1932, -0.7717], \"3\": [1.0333, 0.5966, -0.7717], \"4\": [-1.0333, 0.5966, -0.7717]}, \"atom\": [\"P\", \"H\", \"H\", \"H\"], \"Spin_dn\": 9, \"Spin_up\": 9}',NULL),(52,'CH3SH','CH3SH','{\"xyz\": {\"1\": [-0.85, -0.0344, -0.2], \"2\": [0.9, -0.5125, -0.1219], \"3\": [1.4219, 0.5781, 0.425], \"4\": [-0.9406, 0.8688, -0.8219], \"5\": [-1.4219, -0.8688, -0.6469], \"6\": [-1.2031, 0.1656, 0.821]}, \"atom\": [\"C\", \"S\", \"H\", \"H\", \"H\", \"H\"], \"Spin_dn\": 13, \"Spin_up\": 13}',NULL),(53,'S','S','{\"xyz\": {\"1\": [0, 0, 0]}, \"atom\": [\"S\"], \"Spin_dn\": 7, \"Spin_up\": 9}',NULL),(54,'H2S','H2S','{\"xyz\": {\"1\": [0, 0, 0.103], \"2\": [0, 0.9616, -0.8239], \"3\": [0, -0.9616, -0.8239]}, \"atom\": [\"S\", \"H\", \"H\"], \"Spin_dn\": 9, \"Spin_up\": 9}',NULL),(55,'Na2','Na2','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 3.0789]}, \"atom\": [\"Na\", \"Na\"], \"Spin_dn\": 11, \"Spin_up\": 11}',NULL),(56,'Na','Na','{\"xyz\": {\"1\": [0, 0, 0]}, \"atom\": [\"Na\"], \"Spin_dn\": 5, \"Spin_up\": 6}',NULL),(57,'Si2','Si2','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 2.246]}, \"atom\": [\"Si\", \"Si\"], \"Spin_dn\": 13, \"Spin_up\": 15}',NULL),(58,'P2','P2','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 1.8934]}, \"atom\": [\"P\", \"P\"], \"Spin_dn\": 15, \"Spin_up\": 15}',NULL),(59,'P','P','{\"xyz\": {\"1\": [0, 0, 0]}, \"atom\": [\"P\"], \"Spin_dn\": 6, \"Spin_up\": 9}',NULL),(60,'S2','S2','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 1.8892]}, \"atom\": [\"S\", \"S\"], \"Spin_dn\": 15, \"Spin_up\": 17}',NULL),(61,'S','S','{\"xyz\": {\"1\": [0, 0, 0]}, \"atom\": [\"S\"], \"Spin_dn\": 7, \"Spin_up\": 9}',NULL),(62,'NaCl','NaCl','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 2.3608]}, \"atom\": [\"Na\", \"Cl\"], \"Spin_dn\": 14, \"Spin_up\": 14}',NULL),(63,'CH3Cl','CH3Cl','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 1.781], \"3\": [1.0424, 0, -0.3901], \"4\": [-0.5212, 0.9027, -0.3901], \"5\": [-0.5212, -0.9027, -0.3901]}, \"atom\": [\"C\", \"Cl\", \"H\", \"H\", \"H\"], \"Spin_dn\": 13, \"Spin_up\": 13}',NULL),(64,'TiH','TiH','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 1.7847]}, \"atom\": [\"Ti\", \"H\"], \"Spin_dn\": 11, \"Spin_up\": 12}',NULL),(65,'FCl','FCl','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 0, 1.6283]}, \"atom\": [\"F\", \"Cl\"], \"Spin_dn\": 13, \"Spin_up\": 13}',NULL),(66,'SO2','SO2','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 1.2371, 0.7215], \"3\": [0, -1.2371, 0.7215]}, \"atom\": [\"S\", \"O\", \"O\"], \"Spin_dn\": 16, \"Spin_up\": 16}',NULL),(67,'AlCl3\n','AlCl3','{\"xyz\": {\"1\": [0, 0, 0], \"2\": [0, 2.06, 0.0], \"3\": [1.784, -1.03, 0.0], \"4\": [-1.784, -1.03, 0.0]}, \"atom\": [\"Al\", \"Cl\", \"Cl\", \"Cl\"], \"Spin_dn\": 32, \"Spin_up\": 32}',NULL),(68,'Al','Al','{\"xyz\": {\"1\": [0, 0, 0]}, \"atom\": [\"Al\"], \"Spin_dn\": 6, \"Spin_up\": 7}',NULL);
/*!40000 ALTER TABLE `xyz_molecule` ENABLE KEYS */;
UNLOCK TABLES;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2024-02-06 14:55:40