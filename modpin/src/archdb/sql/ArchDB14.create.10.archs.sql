DROP TABLE IF EXISTS `loop2chain` ;
DROP TABLE IF EXISTS `superloop2chain` ;
DROP TABLE IF EXISTS `loop2superloop` ;
DROP TABLE IF EXISTS `superloop_description` ;
DROP TABLE IF EXISTS `loop_description` ;

CREATE  TABLE IF NOT EXISTS `loop_description` (
  `nid` INT UNSIGNED NOT NULL AUTO_INCREMENT ,
  `loop_id` VARCHAR(15) NOT NULL COMMENT 'composed by the concatenation of pdb, chain and start.' ,
  `type` CHAR(2) NOT NULL ,
  `length` SMALLINT UNSIGNED NOT NULL ,
  `ss1L` TINYINT UNSIGNED NOT NULL COMMENT 'Length of the previous secondary structure' ,
  `ss1moiL` TINYINT UNSIGNED NOT NULL COMMENT 'Length of the previous secondary structure according to its moment of inertia (moi)' ,
  `ss2L` TINYINT UNSIGNED NOT NULL COMMENT 'Length of the following secondary structure' ,
  `ss2moiL` TINYINT UNSIGNED NOT NULL COMMENT 'Length of the following secondary structure according to its moment of inertia (moi)' ,
  `distance` FLOAT(11,6) NOT NULL ,
  `theta` FLOAT(11,6) NOT NULL ,
  `rho` FLOAT(11,6) NOT NULL ,
  `delta` FLOAT(11,6) NOT NULL ,
  `sequence` VARCHAR(255) NOT NULL ,
  `ss` VARCHAR(255) NOT NULL ,
  `exposition` VARCHAR(255) NOT NULL ,
  PRIMARY KEY (`nid`) ,
  UNIQUE INDEX `loop_id_UNIQUE` (`loop_id` ASC) ,
  UNIQUE INDEX `nid_UNIQUE` (`nid` ASC) ,
  INDEX `length_idx` (`length` ASC) ,
  INDEX `type_length_idx` (`type` ASC, `length` ASC) ,
  INDEX `param_idx` (`type` ASC, `length` ASC, `distance` ASC, `theta` ASC, `rho` ASC, `delta` ASC) ,
  INDEX `param_nolength_idx` (`type` ASC, `distance` ASC, `theta` ASC, `rho` ASC, `delta` ASC) )
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `loop2chain` (
  `loop_id` INT UNSIGNED NOT NULL ,
  `chain` INT UNSIGNED NOT NULL ,
  `start` SMALLINT NOT NULL ,
  `idxs` CHAR(1) NULL DEFAULT ' ' ,
  `end` SMALLINT NOT NULL ,
  `idxe` CHAR(1) NULL DEFAULT ' ' ,
  `assignation` ENUM('D','I','H','J') NOT NULL DEFAULT 'D' COMMENT 'D: direct\nI: identical\nH: homology' ,
  `position` SMALLINT UNSIGNED NULL DEFAULT NULL ,
  PRIMARY KEY (`loop_id`, `chain`, `start`) ,
  INDEX `fk_loop2chain_chain2_idx` (`chain` ASC) ,
  INDEX `search_assign_chain` (`assignation` ASC, `chain` ASC, `loop_id` ASC) ,
  INDEX `search_assign_loop` (`assignation` ASC, `loop_id` ASC, `chain` ASC) ,
  CONSTRAINT `fk_loop2chain_loop1`
    FOREIGN KEY (`loop_id` )
    REFERENCES `loop_description` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_loop2chain_chain2`
    FOREIGN KEY (`chain` )
    REFERENCES `chain` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `superloop_description` (
  `nid` INT UNSIGNED NOT NULL AUTO_INCREMENT ,
  `loop_id` VARCHAR(15) NOT NULL COMMENT 'composed by the concatenation of pdb, chain and start.' ,
  `type` CHAR(2) NOT NULL ,
  `length` SMALLINT UNSIGNED NOT NULL ,
  `ss1L` TINYINT UNSIGNED NOT NULL COMMENT 'Length of the previous secondary structure' ,
  `ss1moiL` TINYINT UNSIGNED NOT NULL COMMENT 'Length of the previous secondary structure according to its moment of inertia (moi)' ,
  `ss2L` TINYINT UNSIGNED NOT NULL COMMENT 'Length of the following secondary structure' ,
  `ss2moiL` TINYINT UNSIGNED NOT NULL COMMENT 'Length of the following secondary structure according to its moment of inertia (moi)' ,
  `distance` FLOAT(11,6) NOT NULL ,
  `theta` FLOAT(11,6) NOT NULL ,
  `rho` FLOAT(11,6) NOT NULL ,
  `delta` FLOAT(11,6) NOT NULL ,
  `sequence` VARCHAR(255) NOT NULL ,
  `ss` VARCHAR(255) NOT NULL ,
  `exposition` VARCHAR(255) NOT NULL ,
  `internal_ss_count` TINYINT UNSIGNED NOT NULL DEFAULT 0 ,
  `internal_ss_desc` VARCHAR(45) NULL DEFAULT NULL ,
  PRIMARY KEY (`nid`) ,
  UNIQUE INDEX `loop_id_UNIQUE` (`loop_id` ASC) ,
  UNIQUE INDEX `nid_UNIQUE` (`nid` ASC) ,
  INDEX `length_idx` (`length` ASC) ,
  INDEX `type_length_idx` (`type` ASC, `length` ASC) ,
  INDEX `iss_param_idx` (`internal_ss_count` ASC, `type` ASC, `length` ASC, `distance` ASC, `theta` ASC, `rho` ASC, `delta` ASC) ,
  INDEX `param_idx` (`type` ASC, `length` ASC, `distance` ASC, `theta` ASC, `rho` ASC, `delta` ASC) ,
  INDEX `iss_param_nolength_idx` (`internal_ss_count` ASC, `type` ASC, `distance` ASC, `theta` ASC, `rho` ASC, `delta` ASC) ,
  INDEX `param_nolength_idx` (`type` ASC, `distance` ASC, `theta` ASC, `rho` ASC, `delta` ASC) )
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `superloop2chain` (
  `superloop_id` INT UNSIGNED NOT NULL ,
  `chain` INT UNSIGNED NOT NULL ,
  `start` SMALLINT NULL ,
  `idxs` CHAR(1) NULL ,
  `end` SMALLINT NULL ,
  `idxe` CHAR(1) NULL ,
  `assignation` ENUM('D','I','H') NOT NULL DEFAULT 'D' COMMENT 'D: direct\nI: identical\nH: homology' ,
  `position` SMALLINT UNSIGNED NULL DEFAULT NULL ,
  PRIMARY KEY (`superloop_id`, `chain`,`start`) ,
  INDEX `fk_superloop2chain_chain1_idx` (`chain` ASC) ,
  CONSTRAINT `fk_superloop2chain_superloop_description1`
    FOREIGN KEY (`superloop_id` )
    REFERENCES `superloop_description` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_superloop2chain_chain1`
    FOREIGN KEY (`chain` )
    REFERENCES `chain` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `loop2superloop` (
  `superloop` INT UNSIGNED NOT NULL ,
  `loop` INT UNSIGNED NOT NULL ,
  PRIMARY KEY (`superloop`, `loop`) ,
  INDEX `fk_loop2superloop_loop_description1_idx` (`loop` ASC) ,
  CONSTRAINT `fk_loop2superloop_superloop_description1`
    FOREIGN KEY (`superloop` )
    REFERENCES `superloop_description` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_loop2superloop_loop_description1`
    FOREIGN KEY (`loop` )
    REFERENCES `loop_description` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

DROP procedure IF EXISTS `new_loop2chain`;
DROP procedure IF EXISTS `new_superloop2chain`;

DELIMITER $$

CREATE PROCEDURE new_loop2chain(lid INT UNSIGNED, strnid INT UNSIGNED,
                                lstart SMALLINT, idxs CHAR(1), lend SMALLINT, idxe CHAR(1), 
                                kindof ENUM('D','I','H','J'), lorder SMALLINT UNSIGNED) MODIFIES SQL DATA
BEGIN

  DECLARE finished   INTEGER DEFAULT 0;
  DECLARE correction INTEGER DEFAULT 0;
  DECLARE repechain  INTEGER DEFAULT 0;
  DECLARE list       VARCHAR(255) DEFAULT "";
  DECLARE curchain   CURSOR FOR SELECT r2.chain FROM chainrepeat r1, chainrepeat r2 
                     WHERE r1.chain=strnid AND r1.chain!=r2.chain AND r2.groupid=r1.groupid;

  DECLARE CONTINUE HANDLER FOR NOT FOUND SET finished = 1;
  OPEN curchain;
    IF kindof='D' THEN
      get_ident: LOOP
        FETCH curchain INTO repechain;
        IF finished THEN
          LEAVE get_ident;
        END IF;
        SELECT (c2.start-c1.start) INTO correction 
          FROM chain c1, chain c2 
          WHERE c1.nid=strnid AND c2.nid=repechain;
        INSERT INTO loop2chain(loop_id,chain,start,idxs,end,idxe,assignation) 
          VALUES(lid,repechain,lstart + correction,idxs,lend + correction,idxe,'I');
      END LOOP get_ident;
    END IF;
  IF kindof='H' THEN
      get_ident: LOOP
        FETCH curchain INTO repechain;
        IF finished THEN
          LEAVE get_ident;
        END IF;
        SELECT (c2.start-c1.start) INTO correction 
          FROM chain c1, chain c2 
          WHERE c1.nid=strnid AND c2.nid=repechain;
        INSERT INTO loop2chain(loop_id,chain,start,idxs,end,idxe,assignation) 
          VALUES(lid,repechain,lstart + correction,idxs,lend + correction,idxe,'J');
      END LOOP get_ident;
    END IF;
  CLOSE curchain;

  IF kindof='D' THEN
    INSERT INTO loop2chain(loop_id,chain,start,idxs,end,idxe,assignation,position) 
    VALUES (lid,strnid,lstart,idxs,lend,idxe,kindof,lorder);
  ELSE
    INSERT INTO loop2chain(loop_id,chain,start,idxs,end,idxe,assignation) 
    VALUES (lid,strnid,lstart,idxs,lend,idxe,kindof);
  END IF;
END$$

CREATE PROCEDURE new_superloop2chain(lid INT UNSIGNED, strnid INT UNSIGNED, 
                                     lstart SMALLINT, idxs CHAR(1), lend SMALLINT, idxe CHAR(1), 
                                    kindof ENUM('D','I','H'), lorder SMALLINT UNSIGNED) MODIFIES SQL DATA
BEGIN

  DECLARE finished   INTEGER DEFAULT 0;
  DECLARE correction INTEGER DEFAULT 0;
  DECLARE repechain  INTEGER DEFAULT 0;
  DECLARE list       VARCHAR(255) DEFAULT "";

  DECLARE curchain CURSOR FOR SELECT r2.chain FROM chainrepeat r1, chainrepeat r2 
                   WHERE r1.chain=strnid AND r1.chain!=r2.chain AND r2.groupid=r1.groupid;

  DECLARE CONTINUE HANDLER FOR NOT FOUND SET finished = 1;
  OPEN curchain;

      get_ident: LOOP
        FETCH curchain INTO repechain;
        IF finished THEN
          LEAVE get_ident;
        END IF;
        SELECT c2.start-c1.start INTO correction 
          FROM chain c1, chain c2 
          WHERE c1.nid=strnid AND c2.nid=repechain;
        INSERT INTO superloop2chain(superloop_id,chain,start,idxs,end,idxe,assignation) 
          VALUES(lid,repechain,(lstart + correction),idxs,(lend + correction),idxe,'I');
      END LOOP get_ident;

  CLOSE curchain;

  IF kindof='D' THEN
    INSERT INTO superloop2chain(superloop_id,chain,start,idxs,end,idxe,assignation,position) 
    VALUES (lid,strnid,lstart,idxs,lend,idxe,kindof,lorder);
  ELSE
    INSERT INTO superloop2chain(superloop_id,chain,start,idxs,end,idxe,assignation) 
    VALUES (lid,strnid,lstart,idxs,lend,idxe,kindof);
  END IF;
END$$

DELIMITER ;
