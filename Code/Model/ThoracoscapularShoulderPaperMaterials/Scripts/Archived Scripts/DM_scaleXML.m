import org.opensim.modeling.*
file = uigetfile('*.xml')
%Select scale tool and apply it to the xml file
SCTool = ScaleTool(fullfile(file))
%Getting to measures under model scalar which is under scale tool
ms = SCTool.getModelScaler().getMeasurementSet()
%thorax x bone length
a = ms.get('Thorax_1').getMarkerPairSet()
% a.get('mptest').getPropertyByIndex(0) gives you name of marker
prop = a.get('mptest').updPropertyByIndex(0)
PropertyHelper.setValueString('ai', prop, 0)
a.get('mptest').getPropertyByIndex(0) 
PropertyHelper.setValueString('ts', prop, 1)
a.get('mptest').getPropertyByIndex(0)
b = ms.get('thorax_x').getBodyScaleSet()
b.get('thorax').getPropertyByIndex(0)
prop = b.get('thorax').updPropertyByIndex(0)
PropertyHelper.setValueString('Y', prop, 0)
b.get('thorax').getPropertyByIndex(0)
SCTool.print('test.xml')