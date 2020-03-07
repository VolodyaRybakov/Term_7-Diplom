import json

def getData(station, satelite):
    with open('data.json', 'r', encoding='utf-8') as fh:
        data = json.load(fh)

    a = 0
    try:
        a = (data['satelites'][satelite]['perigee'])/(1-data['satelites'][satelite]['eccentricity'])
    except:
        a = data['other_constants']['earth'] + data['satelites'][satelite]['altitude']

    result = {
        'radial_distance': 6371302 + data['stations'][station]['elevation'],
        'latitude': data['stations'][station]['latitude'],
        'longitude': data['stations'][station]['longitude'],
        'inclination': data['satelites'][satelite]['inclination'],
        'eccentricity': data['satelites'][satelite]['eccentricity'],
        'a': a
    }

    return result

if __name__ == "__main__":
    getData('mendeleevo2', 'cryosat2')