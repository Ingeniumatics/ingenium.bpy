import bpy
import sys
import os


# --- 모듈 파일이 있는 폴더 경로 설정 ---
# 예시: C:/my_blender_scripts 또는 /Users/username/blender_modules
# 스크립트 파일과 같은 폴더에 있다면 아래처럼 상대 경로 사용 가능
#module_dir = os.path.dirname(os.path.abspath(__file__)) # 현재 스크립트 파일의 폴더 경로
# 만약 다른 특정 폴더에 있다면 직접 경로 지정:
module_dir = "D:/ingenium.bpy/src/ingenium/"

if module_dir not in sys.path:
    sys.path.append(module_dir)

# --- 모듈 임포트 ---
try:
    import involute_gear # 파일 이름 (확장자 .py 제외)
    # 개발 중이라면 리로드 기능 추가 (선택 사항)
    import importlib
    importlib.reload(involute_gear)

except ImportError:
    print(f"오류: '{module_dir}' 경로에서 'gear_generator.py' 모듈을 찾을 수 없습니다.")
    # raise # 필요시 에러 발생시켜 중단

# --- 함수 호출 및 사용 ---
if 'involute_gear' in sys.modules:
    gear = involute_gear.InvoluteGear(module=2.0, num_teeth=15, pressure_angle_deg=20.0)
    print("-" * 20)

    # 2. Generate vertices for the full 2D profile
    num_profile_points = 20
    num_fillet_points = 10
    gear_vertices = gear.generate_full_gear_vertices(
        num_profile_points=num_profile_points,
        num_fillet_points=num_fillet_points
    )
