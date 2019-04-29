% Cache knowledge base from database.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 11/26/2012
function [knowledgeBase, knowledgeBaseWID] = cacheKnowledgeBase();

% connect to database
dbConnectionParameters = config();
database = edu.jiangnan.fmme.db.MySQLDatabase(dbConnectionParameters);

% construct latest knowledge base from database
knowledgeBaseWID = edu.jiangnan.fmme.cell.kb.KnowledgeBaseUtil.selectLatestKnowledgeBase(database);
knowledgeBase = edu.jiangnan.fmme.cell.kb.KnowledgeBase(database, knowledgeBaseWID.WID);

%serialize
knowledgeBase.serializeLinks();
knowledgeBaseWID = 1;
% save
save('data/knowledgeBase.mat', 'knowledgeBase', 'knowledgeBaseWID');

%deserialize
knowledgeBase.deserializeLinks();

% cleanup
database.close();