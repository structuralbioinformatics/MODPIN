DROP TABLE IF EXISTS `web_wellcome_summary` ;
DROP TABLE IF EXISTS `web_initial_method_distribution_summary` ;
DROP TABLE IF EXISTS `ws_mainDB` ;

CREATE  TABLE IF NOT EXISTS `web_wellcome_summary` (
  `lcount` BIGINT NOT NULL ,
  `type` CHAR(2) NOT NULL ,
  `cluster` TINYINT(1) NOT NULL ,
  `description` VARCHAR(45) NOT NULL )
ENGINE = InnoDB;

CREATE  TABLE IF NOT EXISTS `web_initial_method_distribution_summary` (
  `count` BIGINT NOT NULL ,
  `type` CHAR(2) NOT NULL ,
  `name` LONGTEXT NULL )
ENGINE = InnoDB;


CREATE  TABLE IF NOT EXISTS `ws_mainDB` (
  `code` INT NOT NULL ,
  `method` VARCHAR(45) NULL ,
  `type` CHAR(2) NOT NULL ,
  `length` VARCHAR(255) NULL ,
  `loop_num` BIGINT NULL ,
  INDEX `fk_ws_mainDB_method1_idx` (`code` ASC) ,
  INDEX `fk_ws_mainDB_type_description1_idx` (`type` ASC) ,
  CONSTRAINT `fk_ws_mainDB_method1`
    FOREIGN KEY (`code` )
    REFERENCES `method` (`nid` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT `fk_ws_mainDB_type_description1`
    FOREIGN KEY (`type` )
    REFERENCES `type_description` (`type` )
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
ENGINE = InnoDB;

SET autocommit = 0;
START TRANSACTION;
INSERT INTO web_wellcome_summary (SELECT 
    COUNT(DISTINCT(ld.nid)) as lcount, ld.type, IF(lc.cluster_nid IS NULL,0,1) AS cluster, t.description 
    FROM loop_description ld 
         LEFT JOIN loop2cluster lc ON lc.loop_nid=ld.nid 
         JOIN type_description t ON t.type=ld.type 
    GROUP BY type,cluster);
INSERT INTO web_initial_method_distribution_summary (SELECT COUNT(t1.lcount) AS count, t1.type, IF(t1.method IS NULL,'None',t1.method) AS name 
    FROM (SELECT ld.nid as lcount, ld.type, GROUP_CONCAT(DISTINCT(m.name)) AS method 
            FROM loop_description ld 
            LEFT JOIN loop2cluster lc ON lc.loop_nid=ld.nid 
            LEFT JOIN cluster_subclass cs ON cs.nid=lc.cluster_nid 
            LEFT JOIN cluster_class cc ON cc.nid=cs.class_nid
             LEFT JOIN method m ON cc.method=m.nid GROUP BY ld.nid ORDER BY m.nid) AS t1 
    GROUP BY type, name);
INSERT INTO ws_mainDB (SELECT 
  m.nid AS code, m.name AS method, ld.type, ld.length, count(ld.loop_id) AS loop_num 
  FROM loop_description ld 
    JOIN loop2cluster lc ON lc.loop_nid=ld.nid 
    JOIN cluster_subclass cs ON cs.nid=lc.cluster_nid 
    JOIN cluster_class cc ON cc.nid=cs.class_nid 
    JOIN method m on cc.method=m.nid 
  GROUP BY method,type,length);
COMMIT;

