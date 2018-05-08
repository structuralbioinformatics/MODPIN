DROP TABLE IF EXISTS `loop2cluster` ;
DROP TABLE IF EXISTS `cluster_subclass` ;
DROP TABLE IF EXISTS `cluster_class` ;
DROP TABLE IF EXISTS `type_description` ;
DROP TABLE IF EXISTS `method` ;

CREATE  TABLE IF NOT EXISTS `method` (
  `nid` INT NOT NULL AUTO_INCREMENT ,
  `name` VARCHAR(45) NOT NULL ,
  `description` VARCHAR(255) NOT NULL ,
  `seed` TINYINT(1) NULL ,
  PRIMARY KEY (`nid`) ,
  UNIQUE INDEX `name_UNIQUE` (`name` ASC) )
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `type_description` (
  `type` CHAR(2) NOT NULL ,
  `description` VARCHAR(45) NULL ,
  PRIMARY KEY (`type`) ,
  UNIQUE INDEX `type_UNIQUE` (`type` ASC) )
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `cluster_class` (
  `nid` INT UNSIGNED NOT NULL AUTO_INCREMENT ,
  `method` INT NOT NULL ,
  `type` CHAR(2) NOT NULL ,
  `length` VARCHAR(255) NOT NULL ,
  `class` SMALLINT UNSIGNED NOT NULL ,
  `name` VARCHAR(255) NOT NULL ,
  `size` SMALLINT UNSIGNED NOT NULL ,
  `consensus` VARCHAR(255) NOT NULL ,
  PRIMARY KEY (`nid`) ,
  UNIQUE INDEX `nid_UNIQUE` (`nid` ASC) ,
  INDEX `classification_idx` (`method` ASC, `type` ASC, `length` ASC, `class` ASC) ,
  INDEX `fk_cluster_database1_idx` (`method` ASC) ,
  UNIQUE INDEX `UNIQUE_class` (`method` ASC, `type` ASC, `length` ASC, `class` ASC) ,
  UNIQUE INDEX `UNIQUE_class_name` (`method` ASC, `name` ASC) ,
  INDEX `fk_cluster_class_type_description1_idx` (`type` ASC) ,
  INDEX `name_idx` (`name` ASC) ,
  CONSTRAINT `fk_cluster_database10`
    FOREIGN KEY (`method` )
    REFERENCES `method` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_cluster_class_type_description1`
    FOREIGN KEY (`type` )
    REFERENCES `type_description` (`type` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `cluster_subclass` (
  `nid` INT UNSIGNED NOT NULL AUTO_INCREMENT ,
  `class_nid` INT UNSIGNED NOT NULL ,
  `subclass` SMALLINT UNSIGNED NOT NULL ,
  `name` VARCHAR(255) NOT NULL ,
  `size` SMALLINT UNSIGNED NOT NULL ,
  `seq_pat` VARCHAR(510) NOT NULL ,
  `exp_pat` VARCHAR(510) NOT NULL ,
  `ram_pat` VARCHAR(510) NOT NULL ,
  `dist_range` VARCHAR(255) NULL ,
  `delta_range` VARCHAR(255) NULL ,
  `theta_range` VARCHAR(255) NULL ,
  `rho_range` VARCHAR(255) NULL ,
  `dist_range_min` SMALLINT UNSIGNED NULL ,
  `dist_range_max` SMALLINT UNSIGNED NULL ,
  `delta_range_min` SMALLINT UNSIGNED NULL ,
  `delta_range_max` SMALLINT UNSIGNED NULL ,
  `theta_range_min` SMALLINT UNSIGNED NULL ,
  `theta_range_max` SMALLINT UNSIGNED NULL ,
  `rho_range_min` SMALLINT UNSIGNED NULL ,
  `rho_range_max` SMALLINT UNSIGNED NULL ,
  `RMSD` FLOAT(6,3) NULL ,
  PRIMARY KEY (`nid`) ,
  UNIQUE INDEX `nid_UNIQUE` (`nid` ASC) ,
  INDEX `classification_idx` (`subclass` ASC) ,
  INDEX `fk_cluster_subclass_cluster_class1_idx` (`class_nid` ASC) ,
  INDEX `subclass_idx` (`class_nid` ASC, `subclass` ASC) ,
  UNIQUE INDEX `sublcass_UNIQUE` (`class_nid` ASC, `subclass` ASC) ,
  INDEX `name_idx` (`name` ASC) ,
  INDEX `distance_idx` (`dist_range_min` ASC, `dist_range_max` ASC) ,
  INDEX `delta_idx` (`delta_range_min` ASC, `delta_range_max` ASC) ,
  INDEX `theta_idx` (`theta_range_min` ASC, `theta_range_max` ASC) ,
  INDEX `rho_idx` (`rho_range_min` ASC, `rho_range_max` ASC) ,
  INDEX `rmsd_idx` (`RMSD` ASC) ,
  CONSTRAINT `fk_cluster_subclass_cluster_class1`
    FOREIGN KEY (`class_nid` )
    REFERENCES `cluster_class` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `loop2cluster` (
  `cluster_nid` INT UNSIGNED NOT NULL ,
  `loop_nid` INT UNSIGNED NOT NULL ,
  `clust_order` SMALLINT UNSIGNED NOT NULL ,
  `score` FLOAT(6,3) UNSIGNED NOT NULL ,
  `seq_ali` VARCHAR(255) NOT NULL ,
  `ss_ali` VARCHAR(255) NOT NULL ,
  `exp_ali` VARCHAR(255) NOT NULL ,
  `ram_ali` VARCHAR(255) NOT NULL ,
  INDEX `fk_loop2cluster_loop_idx` (`loop_nid` ASC) ,
  INDEX `fk_loop2_cluster_cluster_idx` (`cluster_nid` ASC) ,
  PRIMARY KEY (`cluster_nid`, `loop_nid`) ,
  CONSTRAINT `fk_loop2cluster_loop`
    FOREIGN KEY (`loop_nid` )
    REFERENCES `loop_description` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_loop2cluster_cluster`
    FOREIGN KEY (`cluster_nid` )
    REFERENCES `cluster_subclass` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

SET autocommit = 0;
START TRANSACTION;

INSERT INTO method (name, description, seed) VALUES ('DS', 'Density Search', 1);
INSERT INTO method (name, description, seed) VALUES ('MCL', 'Markov Cluster Algorithm', 0);
INSERT INTO type_description (type, description) VALUES ('GG','helix-3:10 - helix-3:10');
INSERT INTO type_description (type, description) VALUES ('GH','helix-3:10 - alpha-helix');
INSERT INTO type_description (type, description) VALUES ('GE','helix-3:10 - beta-strand');
INSERT INTO type_description (type, description) VALUES ('HG','alpha-helix - helix-3:10');
INSERT INTO type_description (type, description) VALUES ('HH','alpha-helix - alpha-helix');
INSERT INTO type_description (type, description) VALUES ('HE','alpha-helix - beta-strand');
INSERT INTO type_description (type, description) VALUES ('EG','beta-strand - helix-3:10');
INSERT INTO type_description (type, description) VALUES ('EH','beta-strand - alpha-helix');
INSERT INTO type_description (type, description) VALUES ('BN','beta-beta hairpin');
INSERT INTO type_description (type, description) VALUES ('BK','beta-beta link');
COMMIT;