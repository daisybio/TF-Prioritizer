import { ComponentFixture, TestBed } from '@angular/core/testing';

import { StatisticalEvaluationComponent } from './statistical-evaluation.component';

describe('StatisticalEvaluationComponent', () => {
  let component: StatisticalEvaluationComponent;
  let fixture: ComponentFixture<StatisticalEvaluationComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ StatisticalEvaluationComponent ]
    })
    .compileComponents();
  });

  beforeEach(() => {
    fixture = TestBed.createComponent(StatisticalEvaluationComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
